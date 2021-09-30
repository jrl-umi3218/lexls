list(APPEND LEXLS_COMPILE_DEFINITIONS "$<TARGET_PROPERTY:lexls,COMPILE_DEFINITIONS>")
list(APPEND LEXLS_COMPILE_DEFINITIONS "$<TARGET_PROPERTY:lexls,INTERFACE_COMPILE_DEFINITIONS>")
list(APPEND LEXLS_COMPILE_FLAGS "$<TARGET_PROPERTY:lexls,COMPILE_FLAGS>")
list(APPEND LEXLS_COMPILE_FLAGS "$<TARGET_PROPERTY:lexls,INTERFACE_COMPILE_FLAGS>")
list(APPEND LEXLS_COMPILE_FLAGS "$<TARGET_PROPERTY:lexls,COMPILE_OPTIONS>")
list(APPEND LEXLS_COMPILE_FLAGS "$<TARGET_PROPERTY:lexls,INTERFACE_COMPILE_OPTIONS>")
list(APPEND LEXLS_INCLUDE_DIRECTORIES "$<TARGET_PROPERTY:lexls,INCLUDE_DIRECTORIES>")
list(APPEND LEXLS_INCLUDE_DIRECTORIES "$<TARGET_PROPERTY:lexls,INTERFACE_INCLUDE_DIRECTORIES>")
list(APPEND LEXLS_LINK_FLAGS "$<TARGET_PROPERTY:lexls,LINK_FLAGS>")

if(NOT DEFINED BUILD_OCTAVE_BINDINGS_WITH_CLANG)
  set(BUILD_OCTAVE_BINDINGS_WITH_CLANG ON)
endif()

if(BUILD_OCTAVE_BINDINGS AND NOT BUILD_OCTAVE_BINDINGS_WITH_CLANG)
  set(BUILD_OCTAVE_BINDINGS OFF)
endif()

if(BUILD_OCTAVE_BINDINGS)
  if(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    find_package(OpenMP REQUIRED)
  endif()
  find_program(OCTAVE_CONFIG octave-config)
  if(NOT OCTAVE_CONFIG)
    message(FATAL_ERROR "BUILD_OCTAVE_BINDINGS is ON but the octave-config program is not installed")
  endif()
  execute_process(COMMAND ${OCTAVE_CONFIG} -p BINDIR
                  OUTPUT_VARIABLE OCTAVE_BIN_PATHS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${OCTAVE_CONFIG} --m-site-dir
                  OUTPUT_VARIABLE OCTAVE_M_SITE_DIR
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  find_program(OCTAVE
               HINTS ${OCTAVE_BIN_PATHS}
               NAMES octave)
  set(MAKE_OCTAVE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/octave/$<CONFIGURATION>")
  configure_file(make_octave.in.m "${CMAKE_CURRENT_BINARY_DIR}/make_octave.cmake.m")
  file(GENERATE
    OUTPUT "${MAKE_OCTAVE_LOCATION}/make_octave.m"
    INPUT "${CMAKE_CURRENT_BINARY_DIR}/make_octave.cmake.m")
  # FIXME This might only work in CMake >= 3.20
  add_custom_command(OUTPUT "${MAKE_OCTAVE_LOCATION}/lexlsi.mex" "${MAKE_OCTAVE_LOCATION}/lexlse.mex"
    COMMAND ${CMAKE_COMMAND} -E chdir "${MAKE_OCTAVE_LOCATION}" ${OCTAVE} --no-gui --no-history --silent make_octave.m
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/lexlse.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/lexlsi.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/make_octave.in.m" "${CMAKE_CURRENT_SOURCE_DIR}/lexls_common.h")
  add_custom_target(lexls-octave-bindings ALL DEPENDS "${MAKE_OCTAVE_LOCATION}/lexlsi.mex" "${MAKE_OCTAVE_LOCATION}/lexlse.mex")
  install(FILES "${MAKE_OCTAVE_LOCATION}/lexlsi.mex" "${MAKE_OCTAVE_LOCATION}/lexlse.mex" DESTINATION "${OCTAVE_M_SITE_DIR}")
endif()

if(BUILD_MATLAB_BINDINGS)
  cmake_minimum_required(VERSION 3.16.3) #see https://gitlab.kitware.com/cmake/cmake/-/issues/19382 and https://gitlab.kitware.com/cmake/cmake/-/merge_requests/4206
  find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM)
  set(Matlab_INSTALL_DIR  "${Matlab_ROOT_DIR}/toolbox/lexls" CACHE PATH "installation directory for Matlab binding libraries")
  matlab_add_mex(NAME matlab_lexlse
                 SRC "${CMAKE_CURRENT_SOURCE_DIR}/lexlse.cpp"
                 OUTPUT_NAME lexlse
                 LINK_TO lexls)
  matlab_add_mex(NAME matlab_lexlsi
                 SRC "${CMAKE_CURRENT_SOURCE_DIR}/lexlsi.cpp"
                 OUTPUT_NAME lexlsi
                 LINK_TO lexls)
  # FIXME Maybe not the correct location
  install(TARGETS matlab_lexlse matlab_lexlsi DESTINATION ${Matlab_INSTALL_DIR})
endif()

macro(add_binding_test NAME FILE_IN)
  set(FILE "${CMAKE_CURRENT_SOURCE_DIR}/${FILE_IN}")
  if(BUILD_OCTAVE_BINDINGS)
    add_test(NAME octave-${NAME}
             COMMAND ${CMAKE_COMMAND} -E chdir "${MAKE_OCTAVE_LOCATION}" ${OCTAVE} --no-gui --no-history "${FILE}")
  endif()
  if(BUILD_MATLAB_BINDINGS)
    if(DEFINED CMAKE_CONFIGURATION_TYPES)
      matlab_add_unit_test(NAME ${NAME} UNITTEST_FILE "${FILE}" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIGURATION>")
    else()
      matlab_add_unit_test(NAME ${NAME} UNITTEST_FILE "${FILE}" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    endif()
  endif()
endmacro()

add_binding_test(example_lexlse examples/example_lexlse.m)
add_binding_test(example_lexlsi examples/example_lexlsi.m)
