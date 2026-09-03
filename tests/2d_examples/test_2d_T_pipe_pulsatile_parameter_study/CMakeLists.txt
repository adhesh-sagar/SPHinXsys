set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_VERBOSE_MAKEFILE on)

STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
SET(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")

file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_INPUT_PATH})

# Two executables live in this folder, so the sources are listed explicitly
# rather than swept up with aux_source_directory (which would put both main()
# functions into a single target).
#
#   <folder name>            the T-pipe parameter study
#   womersley_channel        the analytical Womersley benchmark used to validate
#                            the scheme, which the T-pipe geometry cannot do (its
#                            inlet profile is prescribed and its inlet channel is
#                            far shorter than the entrance length).
add_executable(${PROJECT_NAME} T_pipe_pulsatile_parameter_study.cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
target_link_libraries(${PROJECT_NAME} sphinxsys_2d)

add_executable(womersley_channel womersley_channel_validation.cpp)
set_target_properties(womersley_channel PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}"
    VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
target_link_libraries(womersley_channel sphinxsys_2d)

# Copy the sweep driver and the Python analysis package next to the binary so that
# a study can be launched straight from the build tree.
#
# This is a separate ALL target rather than a POST_BUILD command on the executable:
# POST_BUILD only fires when the binary relinks, so edits to the Python scripts or
# to run_study.sh would silently fail to redeploy.
add_custom_target(${PROJECT_NAME}_tools ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/run_study.sh ${EXECUTABLE_OUTPUT_PATH}/run_study.sh
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/analysis ${EXECUTABLE_OUTPUT_PATH}/analysis
    COMMENT "Deploying run_study.sh and analysis/ next to the binary")
add_dependencies(${PROJECT_NAME}_tools ${PROJECT_NAME} womersley_channel)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME}
    --Re=100 --alpha=5 --A=0.5 --dp=0.3 --cycles=2 --case=smoke
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
set_target_properties(${PROJECT_NAME} PROPERTIES LABELS "pulsatile, bifurcation, poiseuille, womersley, parameter_study")
