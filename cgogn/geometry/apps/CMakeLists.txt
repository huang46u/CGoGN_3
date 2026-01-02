cmake_minimum_required(VERSION 3.7.2 FATAL_ERROR)

project(cgogn_geometry_examples
	LANGUAGES CXX
)

set(CGOGN_TEST_PREFIX "test_")

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)

find_package(Torch REQUIRED)
if (Torch_FOUND)
	message(STATUS "Torch found")
	add_executable(ndf_skeleton ndf_skeleton.cpp)
	target_link_libraries(ndf_skeleton
		cgogn::core
		cgogn::ui
		cgogn::io
		cgogn::rendering
		${CMAKE_DL_LIBS}
		${TORCH_LIBRARIES}
	)
	if (WIN32)
	  get_filename_component(_TORCH_PREFIX "${Torch_DIR}/../../.." ABSOLUTE)

	  if (EXISTS "${_TORCH_PREFIX}/lib")
		add_custom_command(TARGET ndf_skeleton POST_BUILD
		  COMMAND ${CMAKE_COMMAND} -E copy_directory
			"${_TORCH_PREFIX}/lib"
			"$<TARGET_FILE_DIR:ndf_skeleton>"
		  COMMENT "Copying Torch DLLs from lib/"
		)
	  endif()

	  if (EXISTS "${_TORCH_PREFIX}/bin")
		add_custom_command(TARGET ndf_skeleton POST_BUILD
		  COMMAND ${CMAKE_COMMAND} -E copy_directory
			"${_TORCH_PREFIX}/bin"
			"$<TARGET_FILE_DIR:ndf_skeleton>"
		  COMMENT "Copying Torch DLLs from bin/"
		)
	  endif()
	endif()

else()
	message(FATAL_ERROR "Torch not found")
endif()

add_executable(filtering filtering.cpp)
target_link_libraries(filtering
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
)

add_executable(heat_method heat_method.cpp)
target_link_libraries(heat_method
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
)

add_executable(sphere_fitting sphere_fitting.cpp)
target_link_libraries(sphere_fitting
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	Blossom5
	${CMAKE_DL_LIBS}
)

add_executable(registration registration.cpp)
target_link_libraries(registration
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	simpleICP
	${CMAKE_DL_LIBS}
)

add_executable(constriction constriction.cpp)
target_link_libraries(constriction
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
)

if(APPLE)
	find_library(CORE_FOUNDATION CoreFoundation)
	find_library(CARBON Carbon)
	target_link_libraries(filtering
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(heat_method
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(registration
		${CORE_FOUNDATION}
		${CARBON}
	)
endif()
if (CGOGN_USE_CGAL)
	find_package(CGAL REQUIRED)
	add_executable(udf_training udf_training.cpp)
	target_link_libraries(udf_training
		CGAL::CGAL
		cgogn::core
		cgogn::ui
		cgogn::io
		cgogn::rendering
		${CMAKE_DL_LIBS}
	)
endif()
