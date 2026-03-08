# Hybrid Optimization And Inside Connectivity Implementation Log

## 1. Background And Goal
- Goal: integrate the inside connectivity workflow from `udf_inside_connectivity` into the current branch, and switch `udf_training.h` to a hybrid optimization pipeline built on surface points plus projected inside points.
- Scope: modify only `cgogn/geometry/ui_modules/udf_training.h`; do not port the changes to `udf_training_alpha_inside.h`.
- Priority: stability first. Performance can be slower if needed, but the full training / top-up / connectivity / skeleton pipeline must stay coherent and debuggable.

## 2. Frozen Decisions
1. Default connectivity mode is Mode A: Mixed KNN.
2. Sphere parameter updates always use Mode A mixed-cloud SQEM + line-quadric clustering.
3. Mode B affects only neighborhood / skeleton / split connectivity, not sphere parameter updates.
4. Inside projected points are refreshed after sampling, filtering, or top-up insertions; they are not fully reprojected every iteration.
5. Top-up in the training loop is optional and enabled by default.
6. Spheres that still have zero inside support after top-up are deleted at end of iteration; if all spheres are deleted, training stops with an error.
7. Empty-sphere deletion follows the active clustering mode: Mode A uses mixed clusters, Mode B uses inside ownership.
8. Build Skeleton keeps two modes: current connectivity mode, and explicit inside-connectivity build.
9. Fitting and initialization must operate on the mixed point cloud when available.
10. The implementation log must be kept up to date as changes are made.

## 3. Change Summary
### Data Structures
- Added mixed optimization cloud `opt_mesh_` with position / normal / area / knn / quadric / line_quadric / sphere / error / color / normal_color / medial-axis attributes.
- Added mixed runtime caches: `opt_kdtree_`, `opt_ma_kdtree_`, `opt_winding_number_`, `opt_wn_bvh_`, and related vectors.
- Added mixed dirty flag `mixed_cloud_dirty_`.
- Added local-cluster connectivity refresh parameter `local_cluster_connectivity_refresh_interval_`.
- Added sphere visualization attribute `spheres_topup_zero_inside_color_`.

### Training Loop
- Added `prepare_hybrid_iteration_data` to unify inside ownership refresh, optional top-up, and mixed-cloud readiness.
- Refactored fitting-data access through `get_sphere_fit_data`, which now selects between `opt_mesh_` and `samples_mesh_`.
- Refactored `build_kdtree`, `compute_samples_area`, `compute_winding_numbers`, `compute_quadrics`, `recompute_samples_normals_pca`, and `compute_initial_medial_axis` to support direct mixed-cloud execution.
- `init_spheres_from_samples` now initializes from the current fitting context, which means mixed cloud is used automatically when ready.

### Connectivity And Skeleton
- Kept both connectivity modes: Mode A (Mixed KNN) and Mode B (Inside Power).
- Added periodic connectivity refresh for local / neighbor clustering instead of rebuilding every iteration.
- Auto-split and skeleton building now use the current connectivity mode consistently.

### Top-Up And Inside Workflow
- Brought in the staged top-up strategy and logging style from `udf_inside_connectivity`.
- Added a dedicated zero-start prepass for spheres with no inside support.
- Added model-dependent projected-point filtering for UDF and MF.
- Added a dedicated sphere color attribute for spheres that still have zero inside ownership after top-up.

## 4. Implementation Records
- [Step 1] Created the implementation log and froze the overall hybrid optimization plan.
- [Step 2] Extended `PointsParameters` with mixed-cloud data, training controls, and pending-delete state.
- [Step 3] Implemented `rebuild_mixed_optimization_cloud` to merge surface samples and projected inside points into `opt_mesh_`.
- [Step 4] Switched the fitting-data pipeline to use mixed cloud when available.
- [Step 5] Changed sphere initialization to operate on the active fitting context instead of hard-coded `samples_*` data.
- [Step 6] Refactored the fitting-data support chain (`build_kdtree`, area, winding number, quadrics, normals, medial axis) to support mixed-cloud execution and mixed runtime caches.
- [Step 7] Added top-up source logs to distinguish manual top-up, training-loop top-up, and skeleton-build top-up.
- [Step 8] Added projected-point filtering based on model type to reject points that drift away from the target alpha level set.
- [Step 9] Stabilized mixed-cloud MA flip pruning by suppressing direct point deletion in the mixed path.
- [Step 10] Added staged `update_spheres` logs and exception reporting for background training.
- [Step 11] Fixed local clustering collapse-to-one-sphere by replacing invalid-seed assignment-to-sphere-0 with per-point global fallback search.
- [Step 12] Added `mixed_cloud_dirty_` and changed the pipeline to rebuild mixed cloud only when surface or inside data actually changes.
- [Step 13] Added normal-color and cluster-color visualization support for `opt_mesh_`.
- [Step 14] Aligned top-up strategy with `udf_inside_connectivity`: cluster bbox build, zero-start prepass, multi-pass fair-budget sampling, and stalled-pass exit.
- [Step 15] Reduced unnecessary per-iteration inside reassignment and switched end-of-iteration deletion to mode-dependent empty-cluster cleanup.
- [Step 16] Removed unconditional connectivity rebuild at the end of every training iteration.
- [Step 17] Added periodic connectivity refresh for neighbor mode through `local_cluster_connectivity_refresh_interval_`.
- [Step 18] Added detailed zero-start prepass diagnostics: per-sphere attempts, accepts, reject categories, and bbox metrics.
- [Step 19] Added `topup_zero_inside_color` to mark spheres that still have zero inside ownership after top-up.
- [Step 20] Fixed top-up log semantics by splitting the printed identifier into `sphere_rank` and `sphere_mesh_index`.
- [Step 21] Added `prepass_reached` and `skipped_due_budget` to zero-start logs so it is clear whether a sphere was actually sampled or simply never reached because the global prepass budget was exhausted.
- [Step 22] Added `cluster_points_raw`, `cluster_points_effective`, and `bbox_source` so zero-start logs distinguish true cluster support from forced fallback activation.
- [Step 23] Added `points_provider_->set_mesh_bb_vertex_position(*p.opt_mesh_, p.opt_position_)` after mixed-cloud rebuild so provider state stays in sync with mixed-cloud geometry, visibility, and point count reporting.
- [Step 24] Rebuilt this log file after a text-encoding failure corrupted non-ASCII content written from the shell path.

## 5. Validation Notes
### Functional Checks
- Mixed cloud is now used by fitting-data computation, sphere initialization, clustering, sphere updates, and error evaluation.
- Mode A and Mode B can both be selected, and skeleton building follows the chosen mode.
- Top-up supports manual triggering, training-loop triggering, and skeleton-build triggering, all with source logs.
- Spheres that remain with zero inside ownership now have a dedicated visualization attribute.

### Performance And Stability Notes
- Unconditional per-iteration connectivity rebuild was identified as a major extra cost and removed from the default loop.
- Neighbor-mode connectivity is now refreshed periodically instead of every iteration.
- Local clustering can still degrade badly if mixed-cloud point-level sphere seeds are invalid and fallback becomes global for all points.
- Mixed-cloud dirty-flagging removed repeated full recomputation of KNN / area / quadrics / medial-axis data when the underlying point cloud is unchanged.

### Diagnostics
- Top-up source is explicitly logged.
- Zero-start prepass now distinguishes failure causes from budget exhaustion.
- Sphere-update stages are logged clearly enough to isolate crashes or stalls.

## 6. Open Issues And Follow-Up
1. `compute_clusters_local` still needs a better mixed-cloud sphere-seed strategy so it does not fall back to global search for all points.
2. Zero-start top-up still relies on rejection sampling and may fail badly for tiny power cells or aggressive Poisson constraints.
3. Mixed-cloud sphere initialization may create far more spheres than the older pure-surface version because the coverage rule is still very local and KNN-driven.
4. Any future automated rewrite of this log file must stay ASCII-only on the current shell path, or use a confirmed UTF-8-safe write path.
