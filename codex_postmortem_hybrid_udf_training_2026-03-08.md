# Codex Postmortem: Hybrid UDF Training Integration
Date: 2026-03-08
Scope: `cgogn/geometry/ui_modules/udf_training.h` and related workflow/logging changes

## 1) Confirmed Mistakes So Far
1. Requirement misunderstanding on projection bbox:
I changed projection bbox handling in a way you did not ask for, then had to revert.
Root cause: interpreted "0-1 issue" as "force 0-1" instead of "fix model-dependent bbox logic".
Prevention: for model-dependent behavior (UDF vs MF), first restate expected mapping explicitly before coding.

2. Temporary context binding hack (`bind_mixed_as_fitting_context`):
I used temporary rebinding to reuse sample-based pipeline instead of refactoring functions to operate directly on mixed cloud.
Root cause: short-term patch over architecture correction.
Prevention: no temporary alias/rebinding for core data path; refactor call chain owner explicitly (`samples_*` vs `opt_*`).

3. Incomplete migration of fitting pipeline ownership:
`build_kdtree / area / winding / quadrics / initial medial axis` remained implicitly sample-bound at one stage.
Root cause: migration done incrementally without an ownership checklist.
Prevention: enforce one checklist for every stage function and verify all data sources are the intended mesh.

4. Top-up loop switch was not strictly respected:
Top-up related path still executed in cases where "top-up in optimization" was expected off.
Root cause: gating condition not centralized.
Prevention: single gate function for top-up execution; all entry points must check same boolean policy.

5. Wrong deletion criterion in some iterations:
Used inside-count based deletion where final agreed rule is mode-dependent cluster-empty deletion.
Root cause: old policy logic remained in loop.
Prevention: one deletion policy function with explicit mode A/mode B rule and no duplicated branches.

6. Compile break (`marked_for_delete` undeclared):
Refactor introduced symbol mismatch.
Root cause: variable rename without full compile-time reference sweep.
Prevention: after each refactor chunk, run symbol grep and immediate local compile (or at least static scan).

7. Severe performance regression in update loop:
Per-iteration rebuild/reassign/connectivity refresh was too frequent and expensive compared with expected behavior.
Root cause: correctness-first patching without cadence controls.
Prevention: cache invalidation policy:
only rebuild mixed cloud when inside points actually changed;
connectivity update every N iterations (configurable).

8. Diagnostics initially insufficient for zero-start prepass failure:
Logs lacked enough per-sphere rejection reasons and bbox validity signals.
Root cause: coarse aggregate logging only.
Prevention: keep per-sphere debug counters and explicit reject taxonomy (`input/owner/poisson/budget/bbox`).

9. Sphere index/log identity confusion:
Logs showed sphere ids far larger than expected sphere count, causing debugging confusion.
Root cause: mixed use of rank/index/mesh-id without clear label.
Prevention: always print all three when debugging: rank, mesh index, and stable id (if any).

10. Initialization stopping policy mismatch:
Init sphere could continue toward requested max even when no uncovered points remained.
Root cause: no hard stop condition tied to "coverable points exhausted".
Prevention: terminate init early when candidate coverage falls below threshold / no new coverage gain.

11. Provider visibility gap (`opt_mesh`):
`opt_mesh` state/count not visible in provider at one stage, reducing observability.
Root cause: data update not mirrored to provider BB/attributes consistently.
Prevention: provider sync checklist for every new mesh: bb, attrs, emit/update hooks.

12. Log file encoding corruption:
Chinese content became `?` after shell-based write path.
Root cause: unsafe write method for UTF-8 text.
Prevention: avoid shell overwrite for non-ASCII logs; use UTF-8-safe patch/write path only.

## 2) Process Failures Behind These Mistakes
1. I optimized for speed of patching over invariant control in several steps.
2. I changed behavior before locking exact trigger conditions (especially top-up and rebuild cadence).
3. I under-specified observability early, which delayed root-cause isolation.

## 3) Hard Rules To Prevent Repeat
1. No behavioral change before a written invariant list exists for that stage.
2. No duplicated policy logic for top-up/delete/mode switch; centralize each into one function.
3. No temporary data rebinding for core pipeline ownership.
4. No large-loop recomputation unless cache-invalid condition is explicit and logged.
5. Every heavy stage must log start/end + cardinalities + trigger reason.
6. Any non-ASCII file write must use UTF-8-safe edit path only.

## 4) Execution Checklist I Will Follow
1. Before coding:
Write 5-line invariants for the stage being changed.
2. During coding:
Keep one authoritative function each for top-up gating, rebuild gating, delete policy.
3. After coding:
Verify with grep that no stale policy branches remain.
4. Before handoff:
Confirm logs include counts for surface/inside/mixed/spheres/deleted and trigger reasons.
5. Before editing docs/logs with Chinese:
Use UTF-8-safe patch path only; never shell overwrite.

## 5) Immediate Commit Discipline For This Branch
1. Separate commits by intent:
policy fix, performance cadence, logging, visualization.
2. Each commit message must include:
what changed, why, and which invariant it enforces.
3. If behavior differs from prior branch:
document as "intentional delta" in implementation log.
