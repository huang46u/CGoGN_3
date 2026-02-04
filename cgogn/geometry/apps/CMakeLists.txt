cmake_minimum_required(VERSION 3.7.2 FATAL_ERROR)

project(cgogn_geometry_examples
	LANGUAGES CXX
)

set(CGOGN_TEST_PREFIX "test_")

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)

#find_package(Torch REQUIRED)
#if (Torch_FOUND)
#	message(STATUS "Torch found")
#	add_executable(ndf_skeleton ndf_skeleton.cpp)
#	target_link_libraries(ndf_skeleton
#		cgogn::core
#		cgogn::ui
#		cgogn::io
#		cgogn::rendering
#		${CMAKE_DL_LIBS}
#		${TORCH_LIBRARIES}
#	)
#	if (WIN32)
#	  get_filename_component(_TORCH_PREFIX "${Torch_DIR}/../../.." ABSOLUTE)
#
#	  if (EXISTS "${_TORCH_PREFIX}/lib")
#		add_custom_command(TARGET ndf_skeleton POST_BUILD
#		  COMMAND ${CMAKE_COMMAND} -E copy_directory
#			"${_TORCH_PREFIX}/lib"
#			"$<TARGET_FILE_DIR:ndf_skeleton>"
#		  COMMENT "Copying Torch DLLs from lib/"
#		)
#	  endif()
#
#	  if (EXISTS "${_TORCH_PREFIX}/bin")
#		add_custom_command(TARGET ndf_skeleton POST_BUILD
#		  COMMAND ${CMAKE_COMMAND} -E copy_directory
#			"${_TORCH_PREFIX}/bin"
#			"$<TARGET_FILE_DIR:ndf_skeleton>"
#		  COMMENT "Copying Torch DLLs from bin/"
#		)
#	  endif()
#	endif()
#
#else()
#	message(FATAL_ERROR "Torch not found")
#endif()
#
#if (CGOGN_USE_CGAL)
#	find_package(CGAL REQUIRED)
#	find_package(Torch REQUIRED)
#	if(Torch_FOUND)
#		message(STATUS "Torch found for udf_training")
#		add_executable(udf_training udf_training.cpp)
#		target_link_libraries(udf_training
#			CGAL::CGAL
#			cgogn::core
#			cgogn::ui
#			cgogn::io
#			cgogn::rendering
#			${CMAKE_DL_LIBS}
#			${TORCH_LIBRARIES}
#		)
#		if (WIN32)
#		  get_filename_component(_TORCH_PREFIX "${Torch_DIR}/../../.." ABSOLUTE)
#
#		  if (EXISTS "${_TORCH_PREFIX}/lib")
#			add_custom_command(TARGET udf_training POST_BUILD
#			  COMMAND ${CMAKE_COMMAND} -E copy_directory
#				"${_TORCH_PREFIX}/lib"
#				"$<TARGET_FILE_DIR:udf_training>"
#			  COMMENT "Copying Torch DLLs from lib/"
#			)
#		  endif()
#
#		  if (EXISTS "${_TORCH_PREFIX}/bin")
#			add_custom_command(TARGET udf_training POST_BUILD
#			  COMMAND ${CMAKE_COMMAND} -E copy_directory
#				"${_TORCH_PREFIX}/bin"
#				"$<TARGET_FILE_DIR:udf_training>"
#			  COMMENT "Copying Torch DLLs from bin/"
#			)
#		  endif()
#		endif()
#	else()
#	message(FATAL_ERROR "Torch not found")
#	endif()
#endif()

# ---- Torch (official-style) -------------------------------------------------
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

if (CGOGN_USE_CGAL)
  find_package(CGAL REQUIRED)

  cgogn_add_torch_executable(udf_training
    SOURCES udf_training.cpp
    EXTRA_LIBS CGAL::CGAL
  )

  cgogn_add_torch_executable(udf_training_alpha_inside
    SOURCES udf_training_alpha_inside.cpp
    EXTRA_LIBS CGAL::CGAL
  )
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

