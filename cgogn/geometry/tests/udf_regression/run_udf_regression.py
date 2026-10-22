"""Regression test for the UDF skeleton reconstruction (udf_benchmark).

Runs every case of cases.json (input paths relative to the repository root, fixed seeds) one at a time and compares
the skeleton with golden.json. Since the sphere clusters are filled in sample order, a run is bit-reproducible and
independent of the thread count, so the expected result is an identical skeleton file (same SHA-256). When the hash
differs, the skeleton is compared with the golden PLY kept in the output directory: same counts and connectivity and
vertex/radius differences <= --tolerance is reported as NUMERIC, anything else as FAIL.

    python run_udf_regression.py            # compare against golden.json
    python run_udf_regression.py --update   # record golden.json and the golden PLYs from the current build

Runs at below-normal priority on --cores (default 4-11) with matching OpenMP/MKL thread counts, to keep the machine
usable; the result does not depend on the core count."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
CREATE_SUSPENDED = 0x00000004
COUNT_KEYS = ("sample_points", "final_spheres", "skeleton_vertices", "skeleton_edges", "skeleton_faces")


def parse_cores(text):
    lo, _, hi = text.partition("-")
    return list(range(int(lo), int(hi or lo) + 1))


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def read_ply(path):
    head, body = open(path, "rb").read().split(b"end_header\n", 1)
    counts = {l.split()[1]: int(l.split()[2]) for l in head.decode().splitlines() if l.startswith("element")}
    lines = body.decode().splitlines()
    nv = counts["vertex"]
    vertices = [[float(x) for x in l.split()] for l in lines[:nv]]
    return vertices, lines[nv:]


def run_case(exe, name, case, out_dir, cores):
    cfg = json.loads(json.dumps(case))
    for key in ("input_path", "neural_udf_model_path"):
        if cfg["input"].get(key):
            cfg["input"][key] = os.path.join(REPO, cfg["input"][key]).replace("\\", "/")
    stem = os.path.join(out_dir, name)
    cfg["benchmark"] = {"num_measure_runs": 1, "verbose": False}
    cfg["output"] = {"skeleton_ply": stem + ".ply", "timing_json": stem + "_timing.json"}
    for f in (stem + ".ply", stem + "_timing.json"):
        if os.path.exists(f):
            os.remove(f)
    with open(stem + "_cfg.json", "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=1)
    env = {k: v for k, v in os.environ.items() if not k.startswith("UDF_")}
    env.update({k: str(len(cores)) for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")})
    flags = 0
    if os.name == "nt":
        flags = subprocess.BELOW_NORMAL_PRIORITY_CLASS | CREATE_SUSPENDED
    start = time.time()
    with open(stem + ".log", "w", encoding="utf-8", errors="replace") as log:
        p = subprocess.Popen([exe, stem + "_cfg.json"], cwd=os.path.dirname(exe), stdout=log, stderr=subprocess.STDOUT,
                             env=env, creationflags=flags)
        try:
            import psutil
            proc = psutil.Process(p.pid)
            proc.cpu_affinity(cores)
            if os.name == "nt":
                proc.resume()
        except ImportError:
            if os.name == "nt":
                raise SystemExit("psutil is required on Windows (the process is started suspended)")
        rc = p.wait()
    elapsed = time.time() - start
    if rc != 0 or not os.path.exists(stem + ".ply"):
        return None, elapsed
    timing = json.load(open(stem + "_timing.json"))
    counts = {k: int(timing["counts"][k]) for k in COUNT_KEYS if k in timing["counts"]}
    return {"sha256": sha256(stem + ".ply"), "counts": counts}, elapsed


def compare_numeric(ply, golden_ply, tolerance):
    va, ca = read_ply(ply)
    vb, cb = read_ply(golden_ply)
    if len(va) != len(vb) or ca != cb:
        return False, "connectivity differs"
    diff = max((abs(x - y) for ra, rb in zip(va, vb) for x, y in zip(ra, rb)), default=0.0)
    return diff <= tolerance, f"max vertex/radius difference {diff:.3e}"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--exe", default=os.path.join(REPO, "build", "stage", "bin", "Release", "udf_benchmark.exe"))
    parser.add_argument("--out", default=os.path.join(REPO, "build", "udf_regression"))
    parser.add_argument("--cores", default="4-11", help="logical cores to run on, e.g. 4-11")
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--update", action="store_true", help="record golden.json and the golden PLYs")
    parser.add_argument("--only", nargs="*", help="run only these cases")
    args = parser.parse_args()

    cases = json.load(open(os.path.join(HERE, "cases.json"), encoding="utf-8"))
    if args.only:
        cases = {k: v for k, v in cases.items() if k in args.only}
    golden_path = os.path.join(HERE, "golden.json")
    golden = json.load(open(golden_path, encoding="utf-8")) if os.path.exists(golden_path) else {}
    run_dir = os.path.join(args.out, "run")
    golden_dir = os.path.join(args.out, "golden")
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(golden_dir, exist_ok=True)
    cores = parse_cores(args.cores)

    failures = 0
    for name, case in cases.items():
        result, elapsed = run_case(args.exe, name, case, run_dir, cores)
        if result is None:
            print(f"FAIL     {name}: udf_benchmark failed, see {os.path.join(run_dir, name)}.log")
            failures += 1
            continue
        if args.update:
            golden[name] = result
            shutil.copyfile(os.path.join(run_dir, name + ".ply"), os.path.join(golden_dir, name + ".ply"))
            print(f"RECORDED {name}: {result['counts']} ({elapsed:.1f} s)")
            continue
        expected = golden.get(name)
        if expected is None:
            print(f"FAIL     {name}: no golden result, run with --update first")
            failures += 1
        elif result["sha256"] == expected["sha256"]:
            print(f"PASS     {name}: identical skeleton ({elapsed:.1f} s)")
        elif result["counts"] != expected["counts"]:
            print(f"FAIL     {name}: counts {result['counts']} != golden {expected['counts']}")
            failures += 1
        else:
            golden_ply = os.path.join(golden_dir, name + ".ply")
            if not os.path.exists(golden_ply):
                print(f"FAIL     {name}: skeleton differs (same counts), golden PLY missing for a numeric comparison")
                failures += 1
                continue
            ok, detail = compare_numeric(os.path.join(run_dir, name + ".ply"), golden_ply, args.tolerance)
            print(f"{'NUMERIC ' if ok else 'FAIL    '} {name}: same counts, {detail}")
            failures += 0 if ok else 1
    if args.update:
        with open(golden_path, "w", encoding="utf-8", newline="\n") as f:
            json.dump(golden, f, indent=2, sort_keys=True)
            f.write("\n")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
