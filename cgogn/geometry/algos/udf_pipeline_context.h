#ifndef CGOGN_GEOMETRY_ALGOS_UDF_PIPELINE_CONTEXT_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_PIPELINE_CONTEXT_H_

#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap0.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/udf_training.h>
#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

namespace cgogn
{

namespace geometry
{

template <typename RaySamplerTag>
struct UDFPipelineContext
{
	using Surface = cgogn::CMap2;
	using Points = cgogn::CMap0;
	using NonManifold = cgogn::IncidenceGraph;
	using Training = cgogn::ui::UDFTraining<Surface, Points, NonManifold, RaySamplerTag>;

	ui::App app;
	ui::MeshProvider<Surface> surface_provider;
	ui::MeshProvider<Points> points_provider;
	ui::MeshProvider<NonManifold> non_manifold_provider;
	ui::SurfaceRender<Surface> surface_render;
	ui::PointCloudRender<Points> point_cloud_render;
	ui::SurfaceRender<NonManifold> non_manifold_render;
	Training udf_training;

	Surface* selected_surface = nullptr;
	Points* selected_points = nullptr;

	UDFPipelineContext()
		: app(true),
		  surface_provider(app),
		  points_provider(app),
		  non_manifold_provider(app),
		  surface_render(app),
		  point_cloud_render(app),
		  non_manifold_render(app),
		  udf_training(app)
	{
		app.init_modules();
		ui::View* view = app.current_view();
		if (view)
		{
			view->link_module(&surface_provider);
			view->link_module(&points_provider);
			view->link_module(&non_manifold_provider);
			view->link_module(&surface_render);
			view->link_module(&point_cloud_render);
			view->link_module(&non_manifold_render);
			view->link_module(&udf_training);
		}
	}
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_PIPELINE_CONTEXT_H_
