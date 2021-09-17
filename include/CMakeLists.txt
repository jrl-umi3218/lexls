#
# Copyright 2013-2021 CNRS-UM LIRMM, CNRS-AIST JRL, INRIA
#

set(HEADERS_INCLUDE_DIR
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>)

set(HEADERS
  ${HEADERS_INCLUDE_DIR}/lexls/cycling.h
  ${HEADERS_INCLUDE_DIR}/lexls/lexlse.h
  ${HEADERS_INCLUDE_DIR}/lexls/lexls.h
  ${HEADERS_INCLUDE_DIR}/lexls/lexlsi.h
  ${HEADERS_INCLUDE_DIR}/lexls/objective.h
  ${HEADERS_INCLUDE_DIR}/lexls/tools.h
  ${HEADERS_INCLUDE_DIR}/lexls/typedefs.h
  ${HEADERS_INCLUDE_DIR}/lexls/utility.h
  ${HEADERS_INCLUDE_DIR}/lexls/workingset.h
)
set(LEXLS_HEADERS ${HEADERS} PARENT_SCOPE)

add_library(${PROJECT_NAME} INTERFACE)
if(COMMAND target_sources)
  target_sources(${PROJECT_NAME} INTERFACE ${HEADERS})
endif()
target_compile_features(${PROJECT_NAME} INTERFACE cxx_auto_type cxx_constexpr)
if(MSVC)
  target_compile_options(${PROJECT_NAME} INTERFACE /source-charset:utf-8)
endif()
if(TARGET Eigen3::Eigen)
  target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
else()
  target_include_directories(${PROJECT_NAME} SYSTEM INTERFACE "${EIGEN3_INCLUDE_DIR}")
endif()
target_include_directories(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)

install(
    TARGETS ${PROJECT_NAME}
    EXPORT  ${PROJECT_NAME}Targets
)
install(DIRECTORY lexls DESTINATION include FILES_MATCHING PATTERN "*.h")


