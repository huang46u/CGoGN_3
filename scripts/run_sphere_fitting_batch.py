"""Run the headless CGoGN sphere-fitting benchmark and validate its CSV output."""

from __future__ import annotations

import argparse
import csv
import ctypes
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True, help="sphere_fitting_batch executable")
    parser.add_argument("--input-dir", type=Path, required=True, help="directory containing OFF meshes")
    parser.add_argument("--output", type=Path, required=True, help="output CSV")
    parser.add_argument("--timeout", type=int, default=10, help="per-configuration soft timeout in seconds")
    parser.add_argument("--max-iterations", type=int, default=1000)
    parser.add_argument("--surface-samples", type=int, default=300000)
    parser.add_argument("--skeleton-resolution", type=int, default=50)
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--hard-timeout", type=int, default=900, help="per-mesh hard timeout in seconds")
    parser.add_argument("--resume", action="store_true", help="reuse complete per-mesh CSV parts")
    return parser.parse_args()


FIELDNAMES = [
    "mesh",
    "target_spheres",
    "metric",
    "use_line_quadric",
    "auto_split",
    "auto_stop",
    "sphere_correction",
    "timeout_seconds",
    "max_iterations",
    "hausdorff_surface_samples",
    "hausdorff_skeleton_resolution",
    "random_seed",
    "actual_spheres",
    "iterations",
    "surface_to_skeleton",
    "skeleton_to_surface",
    "symmetric",
    "surface_sample_count",
    "skeleton_sample_count",
    "preprocess_seconds",
    "sphere_init_seconds",
    "fit_seconds",
    "eval_seconds",
    "total_seconds",
    "status",
]


def disable_windows_error_dialogs() -> None:
    if sys.platform == "win32":
        sem_failcriticalerrors = 0x0001
        sem_nogpfaulterrorbox = 0x0002
        sem_noopenfileerrorbox = 0x8000
        ctypes.windll.kernel32.SetErrorMode(
            sem_failcriticalerrors | sem_nogpfaulterrorbox | sem_noopenfileerrorbox
        )


def failure_rows(mesh_name: str, args: argparse.Namespace, status: str) -> list[dict[str, str]]:
    rows = []
    for target in (50, 100, 200, 250):
        for metric, use_line in (("sqem_euclidean", "0"), ("sqem_line_quadric", "1")):
            row = {field: "0" for field in FIELDNAMES}
            row.update(
                mesh=mesh_name,
                target_spheres=str(target),
                metric=metric,
                use_line_quadric=use_line,
                auto_split="1",
                auto_stop="1",
                sphere_correction="1",
                timeout_seconds=str(args.timeout),
                max_iterations=str(args.max_iterations),
                hausdorff_surface_samples=str(args.surface_samples),
                hausdorff_skeleton_resolution=str(args.skeleton_resolution),
                random_seed=str(args.seed),
                status=status,
            )
            rows.append(row)
    return rows


def read_part(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_part(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def complete_failure_rows(
    mesh_name: str, args: argparse.Namespace, status: str, completed: list[dict[str, str]]
) -> list[dict[str, str]]:
    completed_by_configuration = {
        (row.get("target_spheres"), row.get("metric")): row
        for row in completed
        if row.get("mesh") == mesh_name
    }
    rows = []
    for failure in failure_rows(mesh_name, args, status):
        key = (failure["target_spheres"], failure["metric"])
        rows.append(completed_by_configuration.get(key, failure))
    return rows


def part_matches(rows: list[dict[str, str]], mesh: Path, args: argparse.Namespace) -> bool:
    if len(rows) != 8 or any(row.get("mesh") != mesh.name for row in rows):
        return False
    expected = {
        "timeout_seconds": str(args.timeout),
        "max_iterations": str(args.max_iterations),
        "hausdorff_surface_samples": str(args.surface_samples),
        "hausdorff_skeleton_resolution": str(args.skeleton_resolution),
        "random_seed": str(args.seed),
    }
    return all(all(row.get(field) == value for field, value in expected.items()) for row in rows)


def write_comparison(output: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "mesh",
        "target_spheres",
        "euclidean_actual_spheres",
        "line_actual_spheres",
        "euclidean_surface_to_skeleton",
        "line_surface_to_skeleton",
        "euclidean_skeleton_to_surface",
        "line_skeleton_to_surface",
        "euclidean_symmetric",
        "line_symmetric",
        "euclidean_surface_sample_count",
        "line_surface_sample_count",
        "euclidean_skeleton_sample_count",
        "line_skeleton_sample_count",
        "euclidean_fit_seconds",
        "line_fit_seconds",
        "euclidean_eval_seconds",
        "line_eval_seconds",
        "euclidean_status",
        "line_status",
    ]
    grouped: dict[tuple[str, str], dict[str, dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault((row["mesh"], row["target_spheres"]), {})[row["metric"]] = row
    comparison_rows = []
    for (mesh, target), metrics in sorted(grouped.items(), key=lambda item: (item[0][0], int(item[0][1]))):
        euclidean = metrics.get("sqem_euclidean", {})
        line = metrics.get("sqem_line_quadric", {})
        comparison_rows.append(
            {
                "mesh": mesh,
                "target_spheres": target,
                "euclidean_actual_spheres": euclidean.get("actual_spheres", "0"),
                "line_actual_spheres": line.get("actual_spheres", "0"),
                "euclidean_surface_to_skeleton": euclidean.get("surface_to_skeleton", "0"),
                "line_surface_to_skeleton": line.get("surface_to_skeleton", "0"),
                "euclidean_skeleton_to_surface": euclidean.get("skeleton_to_surface", "0"),
                "line_skeleton_to_surface": line.get("skeleton_to_surface", "0"),
                "euclidean_symmetric": euclidean.get("symmetric", "0"),
                "line_symmetric": line.get("symmetric", "0"),
                "euclidean_surface_sample_count": euclidean.get("surface_sample_count", "0"),
                "line_surface_sample_count": line.get("surface_sample_count", "0"),
                "euclidean_skeleton_sample_count": euclidean.get("skeleton_sample_count", "0"),
                "line_skeleton_sample_count": line.get("skeleton_sample_count", "0"),
                "euclidean_fit_seconds": euclidean.get("fit_seconds", "0"),
                "line_fit_seconds": line.get("fit_seconds", "0"),
                "euclidean_eval_seconds": euclidean.get("eval_seconds", "0"),
                "line_eval_seconds": line.get("eval_seconds", "0"),
                "euclidean_status": euclidean.get("status", "missing"),
                "line_status": line.get("status", "missing"),
            }
        )
    comparison = output.with_name(f"{output.stem}_comparison.csv")
    with comparison.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(comparison_rows)
    markdown = output.with_name(f"{output.stem}_comparison.md")
    with markdown.open("w", encoding="utf-8") as stream:
        stream.write("# Sphere fitting comparison\n\n")
        stream.write(
            "| mesh | target | E spheres | L spheres | E surface→skeleton | L surface→skeleton | "
            "E skeleton→surface | L skeleton→surface | E symmetric | L symmetric | E surface samples | "
            "L surface samples | E skeleton samples | L skeleton samples | E fit s | L fit s | "
            "E eval s | L eval s | E status | L status |\n"
        )
        stream.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|\n")
        for row in comparison_rows:
            stream.write(
                f"| {row['mesh']} | {row['target_spheres']} | {row['euclidean_actual_spheres']} | "
                f"{row['line_actual_spheres']} | {row['euclidean_surface_to_skeleton']} | "
                f"{row['line_surface_to_skeleton']} | {row['euclidean_skeleton_to_surface']} | "
                f"{row['line_skeleton_to_surface']} | {row['euclidean_symmetric']} | {row['line_symmetric']} | "
                f"{row['euclidean_surface_sample_count']} | {row['line_surface_sample_count']} | "
                f"{row['euclidean_skeleton_sample_count']} | {row['line_skeleton_sample_count']} | "
                f"{row['euclidean_fit_seconds']} | {row['line_fit_seconds']} | {row['euclidean_eval_seconds']} | "
                f"{row['line_eval_seconds']} | {row['euclidean_status']} | {row['line_status']} |\n"
            )
    print(f"Wrote comparison tables to {comparison} and {markdown}")


def main() -> int:
    args = parse_args()
    disable_windows_error_dialogs()
    meshes = sorted(args.input_dir.glob("*.off"))
    if not meshes:
        raise SystemExit(f"no OFF files found in {args.input_dir}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    parts_dir = args.output.parent / f"{args.output.stem}_parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, mesh in enumerate(meshes):
        part = parts_dir / f"{index:04d}_{mesh.name}.csv"
        if args.resume and part.exists():
            cached = read_part(part)
            if part_matches(cached, mesh, args):
                rows.extend(cached)
                print(f"Resume {mesh.name}")
                continue
        command = [
            str(args.exe),
            "--mesh",
            str(mesh),
            "--output",
            str(part),
            "--timeout",
            str(args.timeout),
            "--max-iterations",
            str(args.max_iterations),
            "--surface-samples",
            str(args.surface_samples),
            "--skeleton-resolution",
            str(args.skeleton_resolution),
            "--seed",
            str(args.seed),
        ]
        print("Running:", " ".join(command))
        status = None
        try:
            completed = subprocess.run(command, check=False, timeout=args.hard_timeout or None)
            if completed.returncode != 0:
                status = f"process_exit_{completed.returncode}"
        except subprocess.TimeoutExpired:
            status = f"hard_timeout_{args.hard_timeout}s"
        if status is not None:
            completed_rows = read_part(part) if part.exists() else []
            write_part(part, complete_failure_rows(mesh.name, args, status, completed_rows))
        cached = read_part(part) if part.exists() else []
        if len(cached) != 8:
            write_part(part, failure_rows(mesh.name, args, "incomplete_process_output"))
            cached = read_part(part)
        rows.extend(cached)

    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    expected_rows = len(meshes) * 4 * 2
    if len(rows) != expected_rows:
        raise SystemExit(f"expected {expected_rows} rows, got {len(rows)}")
    statuses = {}
    for row in rows:
        statuses[row["status"]] = statuses.get(row["status"], 0) + 1
    print(f"Wrote {len(rows)} rows to {args.output}")
    print("Statuses:", ", ".join(f"{key}={value}" for key, value in sorted(statuses.items())))
    write_comparison(args.output, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
