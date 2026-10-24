cmake_minimum_required(VERSION 3.7.2 FATAL_ERROR)

project(cgogn_geometry_examples
	LANGUAGES CXX
)

set(CGOGN_TEST_PREFIX "test_")

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)

# ---- Torch (optional, official-style) ---------------------------------------
# Torch is searched for only when CGOGN_USE_TORCH is ON. Without it, every UDF
# executable is built without Torch and with CGOGN_UDF_WITH_TORCH=0.
if (CGOGN_USE_TORCH)
  find_package(Torch REQUIRED)
  message(STATUS "Torch found: ${Torch_DIR}")

  # Official PyTorch C++ docs recommend appending TORCH_CXX_FLAGS
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

  # Helper: add an executable that links CGoGN + Torch (+ optional extras),
  # and copies Torch DLLs next to the exe on Windows/MSVC
  function(cgogn_add_torch_executable target)
    # usage:
    #   cgogn_add_torch_executable(name
    #     SOURCES file1.cpp [file2.cpp ...]
    #     EXTRA_LIBS libA libB ...
    #   )
    set(options)
    set(oneValueArgs)
    set(multiValueArgs SOURCES EXTRA_LIBS)
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    add_executable(${target} ${ARG_SOURCES})

    target_link_libraries(${target}
      PRIVATE
        cgogn::geometry
        cgogn::core
        cgogn::ui
        cgogn::io
        cgogn::rendering
        ${CMAKE_DL_LIBS}
        "${TORCH_LIBRARIES}"
        ${ARG_EXTRA_LIBS}
    )

    set_property(TARGET ${target} PROPERTY CXX_STANDARD 17)
    set_property(TARGET ${target} PROPERTY CXX_STANDARD_REQUIRED ON)

    # On Windows with MSVC, copy libtorch DLLs next to the executable
    if (MSVC)
      file(GLOB TORCH_DLLS "${TORCH_INSTALL_PREFIX}/lib/*.dll")
      if (TORCH_DLLS)
        add_custom_command(TARGET ${target}
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different
                  ${TORCH_DLLS}
                  $<TARGET_FILE_DIR:${target}>
          COMMENT "Copying Torch DLLs to output directory"
        )
      endif()
    endif()
  endfunction()

  # Torch-based executables
  cgogn_add_torch_executable(ndf_skeleton
    SOURCES ndf_skeleton.cpp
  )
endif()

# Helper: UDF executables (udf_training, udf_benchmark, udf_residual_prune_test).
# With Torch: linked against Torch, CGOGN_UDF_WITH_TORCH=1 (neural input available).
# Without Torch: plain executable on the same cgogn libraries, CGOGN_UDF_WITH_TORCH=0.
function(cgogn_add_udf_executable target)
  # usage:
  #   cgogn_add_udf_executable(name
  #     SOURCES file1.cpp [file2.cpp ...]
  #     EXTRA_LIBS libA libB ...
  #   )
  set(options)
  set(oneValueArgs)
  set(multiValueArgs SOURCES EXTRA_LIBS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if (CGOGN_USE_TORCH)
    cgogn_add_torch_executable(${target}
      SOURCES ${ARG_SOURCES}
      EXTRA_LIBS ${ARG_EXTRA_LIBS}
    )
    target_compile_definitions(${target} PRIVATE CGOGN_UDF_WITH_TORCH=1)
  else()
    add_executable(${target} ${ARG_SOURCES})
    target_link_libraries(${target}
      PRIVATE
        cgogn::geometry
        cgogn::core
        cgogn::ui
        cgogn::io
        cgogn::rendering
        ${CMAKE_DL_LIBS}
        ${ARG_EXTRA_LIBS}
    )
    set_property(TARGET ${target} PROPERTY CXX_STANDARD 17)
    set_property(TARGET ${target} PROPERTY CXX_STANDARD_REQUIRED ON)
    target_compile_definitions(${target} PRIVATE CGOGN_UDF_WITH_TORCH=0)
  endif()
endfunction()

if (CGOGN_USE_CGAL)
  find_package(CGAL REQUIRED)

  cgogn_add_udf_executable(udf_training
    SOURCES udf_training.cpp
    EXTRA_LIBS CGAL::CGAL
  )

  cgogn_add_udf_executable(udf_benchmark
    SOURCES
      udf_benchmark.cpp
      ../benchmark/benchmark_config.cpp
      ../benchmark/skeleton_benchmark_runner.cpp
    EXTRA_LIBS CGAL::CGAL
  )

  cgogn_add_udf_executable(udf_stage_sequence_test
    SOURCES
      udf_stage_sequence_test.cpp
      ../benchmark/benchmark_config.cpp
      ../benchmark/skeleton_benchmark_runner.cpp
    EXTRA_LIBS CGAL::CGAL
  )

  cgogn_add_udf_executable(udf_residual_prune_test
    SOURCES udf_residual_prune_test.cpp
    EXTRA_LIBS CGAL::CGAL
  )

  # Self-containment check: every UDF header is compiled alone in its own translation unit.
  # Headers listed in UDF_HEADER_CHECK_EXCLUDE are known not to be self-contained yet.
  set(UDF_HEADER_CHECK_EXCLUDE)
  file(GLOB_RECURSE UDF_CHECKED_HEADERS CONFIGURE_DEPENDS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/.."
    "${CMAKE_CURRENT_SOURCE_DIR}/../algos/udf/*.h")
  list(APPEND UDF_CHECKED_HEADERS
    types/ray_level_set_sampler.h
    types/ray_level_set_sampler_traits.h
  )
  if (NOT CGOGN_USE_TORCH)
    # The neural headers need LibTorch; without it they are not part of the check.
    list(FILTER UDF_CHECKED_HEADERS EXCLUDE REGEX "^algos/udf/neural/")
  endif()
  set(UDF_HEADER_CHECK_SOURCES)
  foreach(header ${UDF_CHECKED_HEADERS})
    if (NOT header IN_LIST UDF_HEADER_CHECK_EXCLUDE)
      string(MAKE_C_IDENTIFIER "${header}" header_id)
      set(check_source "${CMAKE_CURRENT_BINARY_DIR}/udf_header_check/${header_id}.cpp")
      file(CONFIGURE OUTPUT "${check_source}" CONTENT "#include <cgogn/geometry/${header}>
")
      list(APPEND UDF_HEADER_CHECK_SOURCES "${check_source}")
    endif()
  endforeach()
  if (CGOGN_USE_TORCH)
    set(UDF_HEADER_CHECK_TORCH_LIBS "${TORCH_LIBRARIES}")
  else()
    set(UDF_HEADER_CHECK_TORCH_LIBS)
  endif()
  add_library(udf_header_check OBJECT ${UDF_HEADER_CHECK_SOURCES})
  target_link_libraries(udf_header_check PRIVATE cgogn::geometry cgogn::core cgogn::io ${UDF_HEADER_CHECK_TORCH_LIBS} CGAL::CGAL)
  target_compile_definitions(udf_header_check PRIVATE CGOGN_UDF_WITH_TORCH=$<BOOL:${CGOGN_USE_TORCH}>)
  set_property(TARGET udf_header_check PROPERTY CXX_STANDARD 17)
  set_property(TARGET udf_header_check PROPERTY CXX_STANDARD_REQUIRED ON)
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

add_executable(fast_winding_number_test fast_winding_number_test.cpp)
target_link_libraries(fast_winding_number_test
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

add_executable(shapenet_obj_to_normalized_ply shapenet_obj_to_normalized_ply.cpp)
target_link_libraries(shapenet_obj_to_normalized_ply
	cgogn::geometry
	cgogn::core
	cgogn::io
	${CMAKE_DL_LIBS}
)
set_property(TARGET shapenet_obj_to_normalized_ply PROPERTY CXX_STANDARD 17)
set_property(TARGET shapenet_obj_to_normalized_ply PROPERTY CXX_STANDARD_REQUIRED ON)

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


