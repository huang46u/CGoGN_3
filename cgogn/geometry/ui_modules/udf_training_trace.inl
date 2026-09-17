// Included inside UDFTraining. The optimizer owns live state; diagnostics use a deep geometry copy.
void headless_copy_reconstruction_state_prepared(POINTS& source, POINTS& destination)
{
    if (&source == &destination)
        throw std::invalid_argument("Snapshot destination must be independent.");
    const PointsParameters& src = points_parameters_.at(&source);
    PointsParameters& dst = points_parameters_.at(&destination);
    if (src.input_mode_ != INPUT_SURFACE_MESH)
        throw std::runtime_error("Snapshot reconstruction currently supports surface-mesh inputs only.");
    dst.input_mode_ = src.input_mode_;
    dst.iteration_count_ = src.iteration_count_;
    dst.fitting_data_computed_ = true;
    std::unordered_map<uint32, PVertex> samples;
    if (nb_cells<PVertex>(*dst.samples_mesh_) == 0)
    {
        foreach_cell(*src.samples_mesh_, [&](PVertex v) {
            const uint32 i = index_of(*src.samples_mesh_, v);
            PVertex w = add_vertex(*dst.samples_mesh_);
            const uint32 j = index_of(*dst.samples_mesh_, w);
            samples.emplace(i, w);
            (*dst.samples_position_)[j] = (*src.samples_position_)[i];
            (*dst.samples_normal_)[j] = (*src.samples_normal_)[i];
            (*dst.samples_area_)[j] = (*src.samples_area_)[i];
            (*dst.samples_quadric_)[j] = (*src.samples_quadric_)[i];
            (*dst.samples_line_quadric_)[j] = (*src.samples_line_quadric_)[i];
            (*dst.samples_ma_position_)[j] = (*src.samples_ma_position_)[i];
            (*dst.samples_ma_radius_)[j] = (*src.samples_ma_radius_)[i];
            (*dst.samples_normal_color_)[j] = (*src.samples_normal_color_)[i];
            return true;
        });
        auto ids = get_or_add_attribute<uint32, PVertex>(*dst.samples_mesh_, "source_id");
        foreach_cell(*src.samples_mesh_, [&](PVertex v) {
            const uint32 i = index_of(*src.samples_mesh_, v);
            const uint32 j = index_of(*dst.samples_mesh_, samples.at(i));
            (*ids)[j] = i;
            auto& neighbors = (*dst.samples_knn_)[j];
            neighbors.clear();
            for (PVertex n : (*src.samples_knn_)[i])
                neighbors.push_back(samples.at(index_of(*src.samples_mesh_, n)));
            return true;
        });
        build_kdtree(dst);
    }
    else
    {
        auto ids = get_attribute<uint32, PVertex>(*dst.samples_mesh_, "source_id");
        foreach_cell(*dst.samples_mesh_, [&](PVertex v) {
            samples.emplace(value<uint32>(*dst.samples_mesh_, ids, v), v);
            return true;
        });
        if (samples.size() != nb_cells<PVertex>(*src.samples_mesh_))
            throw std::runtime_error("Source sample count changed during optimization.");
    }
    clear(*dst.skeleton_);
    dst.skeleton_faces_map_.clear();
    dst.skeleton_tets_.clear();
    invalidate_skeleton_face_score_cache(dst);
    invalidate_topology_stage_snapshot(dst);
    clear(*dst.spheres_);
    std::unordered_map<uint32, PVertex> spheres;
    auto ids = get_or_add_attribute<uint32, PVertex>(*dst.spheres_, "source_id");
    foreach_cell(*src.spheres_, [&](PVertex v) {
        const uint32 i = index_of(*src.spheres_, v);
        PVertex w = add_vertex(*dst.spheres_);
        const uint32 j = index_of(*dst.spheres_, w);
        spheres.emplace(i, w);
        (*ids)[j] = i;
        (*dst.spheres_position_)[j] = (*src.spheres_position_)[i];
        (*dst.spheres_radius_)[j] = (*src.spheres_radius_)[i];
        (*dst.spheres_color_)[j] = (*src.spheres_color_)[i];
        (*dst.spheres_cluster_color_)[j] = (*src.spheres_cluster_color_)[i];
        (*dst.spheres_cluster_area_)[j] = (*src.spheres_cluster_area_)[i];
        (*dst.spheres_error_)[j] = (*src.spheres_error_)[i];
        (*dst.spheres_sqem_lambda_)[j] = (*src.spheres_sqem_lambda_)[i];
        (*dst.spheres_do_not_split_)[j] = true;
        (*dst.spheres_cluster_)[j].clear();
        (*dst.spheres_neighbor_clusters_)[j].clear();
        for (PVertex s : (*src.spheres_cluster_)[i])
            (*dst.spheres_cluster_)[j].push_back(samples.at(index_of(*src.samples_mesh_, s)));
        return true;
    });
    dst.nb_spheres_ = static_cast<uint32>(spheres.size());
    foreach_cell(*src.samples_mesh_, [&](PVertex v) {
        const uint32 i = index_of(*src.samples_mesh_, v);
        const uint32 j = index_of(*dst.samples_mesh_, samples.at(i));
        const PVertex owner = (*src.samples_sphere_)[i];
        (*dst.samples_sphere_)[j] = owner.is_valid()
            ? spheres.at(index_of(*src.spheres_, owner)) : PVertex();
        (*dst.samples_error_)[j] = (*src.samples_error_)[i];
        return true;
    });
}

void headless_export_training_state_prepared(POINTS& points, const std::string& directory, bool include_samples)
{
    PointsParameters& p = points_parameters_.at(&points);
    const std::filesystem::path dir(directory);
    std::filesystem::create_directories(dir);
    io::PointExportAttributeSelection<POINTS> attrs;
    attrs.vertex_attributes.push_back(p.spheres_radius_);
    auto ids = get_attribute<uint32, PVertex>(*p.spheres_, "source_id");
    if (ids)
        attrs.vertex_attributes.push_back(ids);
    attrs.vertex_color_attribute = p.spheres_cluster_color_;
    points_provider_->save_points_ply_to_file(*p.spheres_, p.spheres_position_.get(),
        (dir / "spheres.ply").string(), attrs);
    auto sample_ids = get_attribute<uint32, PVertex>(*p.samples_mesh_, "source_id");
    std::ofstream assignment(dir / "clusters.csv");
    assignment << "sample_id,sphere_id\n";
    foreach_cell(*p.samples_mesh_, [&](PVertex v) {
        const uint32 i = index_of(*p.samples_mesh_, v);
        const PVertex owner = (*p.samples_sphere_)[i];
        const uint32 sphere_index = owner.is_valid() ? index_of(*p.spheres_, owner) : INVALID_INDEX;
        assignment << (sample_ids ? (*sample_ids)[i] : i) << ","
            << (sphere_index == INVALID_INDEX ? INVALID_INDEX : (ids ? (*ids)[sphere_index] : sphere_index)) << "\n";
        return true;
    });
    if (!assignment)
        throw std::runtime_error("Failed to write cluster assignments.");
    if (include_samples)
    {
        if (!export_samples_mesh_cluster_ply(p, (dir / "offset_samples.ply").string()))
            throw std::runtime_error("Failed to export offset samples.");
        io::PointExportAttributeSelection<POINTS> normals;
        normals.vertex_attributes.push_back(p.samples_normal_);
        if (sample_ids)
            normals.vertex_attributes.push_back(sample_ids);
        points_provider_->save_points_ply_to_file(*p.samples_mesh_, p.samples_position_.get(),
            (dir / "offset_samples_normals.ply").string(), normals);
    }
}

std::size_t headless_optimization_state_fingerprint(POINTS& points) const
{
    const PointsParameters& p = points_parameters_.at(&points);
    std::size_t hash = 0;
    auto combine = [&](std::size_t x) { hash ^= x + 0x9e3779b9 + (hash << 6) + (hash >> 2); };
    foreach_cell(*p.spheres_, [&](PVertex v) {
        const uint32 i = index_of(*p.spheres_, v);
        combine(i);
        for (int axis = 0; axis < 3; ++axis)
            combine(std::hash<Scalar>{}((*p.spheres_position_)[i][axis]));
        combine(std::hash<Scalar>{}((*p.spheres_radius_)[i]));
        for (PVertex n : (*p.spheres_neighbor_clusters_)[i])
            combine(index_of(*p.spheres_, n));
        return true;
    });
    foreach_cell(*p.samples_mesh_, [&](PVertex v) {
        const PVertex owner = (*p.samples_sphere_)[index_of(*p.samples_mesh_, v)];
        combine(owner.is_valid() ? index_of(*p.spheres_, owner) : INVALID_INDEX);
        return true;
    });
    return hash;
}

