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
else()
	message(FATAL_ERROR "Torch not found")
endif()
if (MSVC)
  file(GLOB TORCH_DLLS "${TORCH_INSTALL_PREFIX}/lib/*.dll")
  add_custom_command(TARGET ndf_skeleton
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different
                     ${TORCH_DLLS}
                     $<TARGET_FILE_DIR:ndf_skeleton>)
endif (MSVC)

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
