# Split UDF Training into Sample-Only + Alpha-Inside Executables

## Summary
Create a new executable `udf_training_alpha_inside` that performs optimization **only** on alpha_inside, while the existing `udf_training` becomes **sample-only** (all alpha_inside logic removed). Split the current `AlphaSetSphereDebugger` into two dedicated debuggers (samples vs alpha_inside). Produce `plan.md` at repo root with this plan before implementation.

## Plan.md (before coding)
1. Create `d:\Code\CGoGN_3\plan.md` and paste the exact plan content below (this block), so it can be reviewed and approved.
2. After approval, begin implementation.

## High-Level Changes
- **New executable**: `udf_training_alpha_inside` (alpha_inside optimization only).
- **Existing executable**: `udf_training` becomes sample-only; all alpha_inside logic removed.
- **New UI module**: `UDFTrainingAlphaInside` in `cgogn/geometry/ui_modules/udf_training_alpha_inside.h`.
- **Split debugger**:
  - `AlphaSamplesSphereDebugger` (samples only)
  - `AlphaInsideSphereDebugger` (alpha_inside only)
- **CMake**: add new executable entry in `cgogn/geometry/apps/CMakeLists.txt`.

---

## Detailed Steps (Decision-Complete)

### 1) Create new alpha-inside header
**File**: `cgogn/geometry/ui_modules/udf_training_alpha_inside.h`

- Start from current `udf_training.h` as baseline (copy).
- **Remove sample-optimization path**:
  - Remove `SphereFitSource` enum and all branching logic that switches between samples vs alpha_inside.
  - Fit/optimization data should be **fixed** to alpha_inside projected data.
- **Keep sample-mesh UI/logic needed for alpha_inside sampling**:
  - Keep alpha-level set sampling (samples mesh), winding number building, etc.
  - Keep ?Alpha Inside Sampling? UI and functions.
- **Alpha-only optimization**:
  - All clustering/optimization uses alpha_inside (projected) data.
  - Preserve the ?component-aware update? checkbox behavior you added, but scoped to alpha-only.
- **Keep needed rendering**:
  - Maintain alpha_inside point rendering, sphere rendering, skeleton rendering for alpha_inside.
- **Remove any sample-only optimization**:
  - Remove sample-specific attributes like `samples_quadric_` if only used for optimization.
  - Keep only if needed for alpha-level set sampling and winding number.

> Outcome: This header compiles standalone and does not depend on sample optimization code paths.

### 2) Make sample-only version of udf_training
**File**: `cgogn/geometry/ui_modules/udf_training.h`

- Remove **all alpha_inside** fields, functions, and UI:
  - `alpha_inside_*` attributes, kdtree, projection, sampling, medial axis, inside cluster logic, etc.
  - `Alpha Inside Sampling` UI section.
  - Any alpha_inside-only options (component update checkbox, cluster mode for alpha_inside).
- Keep **only**:
  - Input handling, sample mesh sampling, sample-based fitting/optimization, rendering.
- Ensure build still works with `udf_training.cpp`.

> Outcome: `udf_training` is clean sample-only.

### 3) Split the sphere debugger
Create two new debuggers so each executable has only one, no source-switching:

**A. AlphaSamplesSphereDebugger**
- **File**: `cgogn/geometry/ui_modules/alpha_samples_sphere_debugger.h`
- Derived from current `alpha_set_sphere_debugger.h`, but:
  - Remove all alpha_inside attributes and logic.
  - Use only sample attributes (`samples_position_`, `samples_normal_`, etc.).
  - Clustering uses sample positions and sample quadrics.
  - UI has no source switching.

**B. AlphaInsideSphereDebugger**
- **File**: `cgogn/geometry/ui_modules/alpha_inside_sphere_debugger.h` (new or reintroduced)
- Use alpha_inside data only:
  - Cluster positions = `alpha_inside_position_`
  - Projected positions = `alpha_inside_projected_position_`
  - Optimization uses projected quadric/line quadric
  - Keep connected-component logic (projected KNN)
- UI has no source switching.

**Note**: `alpha_set_sphere_debugger.h` will be left unused; optional cleanup after verification.

### 4) Create new executable: udf_training_alpha_inside
**Files**:
- `cgogn/geometry/apps/udf_training_alpha_inside.cpp` (new)
  - Copy structure from `udf_training.cpp`
  - Replace include:
    - `udf_training.h` ? `udf_training_alpha_inside.h`
    - `alpha_set_sphere_debugger.h` ? `alpha_inside_sphere_debugger.h`
  - Instantiate `UDFTrainingAlphaInside` and `AlphaInsideSphereDebugger`
  - Update window title to something like **"UDF Training (Alpha Inside)"** for clarity
- `cgogn/geometry/apps/CMakeLists.txt`
  - Add new target:
    ```
    cgogn_add_torch_executable(udf_training_alpha_inside
      SOURCES udf_training_alpha_inside.cpp
      EXTRA_LIBS CGAL::CGAL
    )
    ```
  - Keep existing `udf_training` target unchanged (still sample-only).

### 5) Update sample executable to use sample debugger
**File**: `cgogn/geometry/apps/udf_training.cpp`
- Replace `alpha_set_sphere_debugger.h` with `alpha_samples_sphere_debugger.h`
- Instantiate `AlphaSamplesSphereDebugger` instead of `AlphaSetSphereDebugger`

### 6) Sanity passes / compile checks (non-mutating)
- Ensure both targets build:
  - `udf_training`
  - `udf_training_alpha_inside`
- Verify no include conflicts.

---

## Public API / Interface Changes
- New UI module class: `UDFTrainingAlphaInside<...>`
- New debugger classes:
  - `AlphaSamplesSphereDebugger<Points>`
  - `AlphaInsideSphereDebugger<Points>`
- New executable: `udf_training_alpha_inside`
- Removal of alpha_inside UI/logic from sample-only `UDFTraining`

---

## Tests / Verification Scenarios
1. **Build** both executables.
2. **Sample-only app (`udf_training`)**:
   - Load mesh/point cloud
   - Sample alpha-level set
   - Fit spheres with samples
   - No alpha_inside UI appears
3. **Alpha-inside app (`udf_training_alpha_inside`)**:
   - Load mesh/point cloud
   - Sample alpha-level set (needed for winding number)
   - Sample alpha_inside
   - Run optimization on alpha_inside
   - Use AlphaInsideSphereDebugger, press O without crash

---

## Assumptions
- `plan.md` should live in repo root `d:\Code\CGoGN_3\plan.md`.
- Sample-only version removes all alpha_inside UI and code.
- Debugger split uses two new modules; `alpha_set_sphere_debugger.h` is left unused (optional cleanup later).
