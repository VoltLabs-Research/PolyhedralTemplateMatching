#
# Copyright (C) Rodolfo Herrera Hernandez. All rights reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(polyhedral-template-matching VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(
    POLYHEDRAL_TEMPLATE_MATCHING_REFRESH_LOCAL_DEPS
    "Refresh sibling CoreToolkit and StructureIdentification local builds before linking"
    ON
)

set(
    POLYHEDRAL_TEMPLATE_MATCHING_LOCAL_CORETOOLKIT_BUILD_DIR
    "${CMAKE_SOURCE_DIR}/../CoreToolkit/build-local/build/Release"
    CACHE PATH "Local CoreToolkit build directory used during development"
)
set(
    POLYHEDRAL_TEMPLATE_MATCHING_LOCAL_STRUCTURE_IDENTIFICATION_BUILD_DIR
    "${CMAKE_SOURCE_DIR}/../StructureIdentification/build-local/build/Release"
    CACHE PATH "Local StructureIdentification build directory used during development"
)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)

find_package(coretoolkit REQUIRED)
find_package(structure-identification REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)

set(POLYHEDRAL_TEMPLATE_MATCHING_BOOST_TARGET "")
if(TARGET Boost::headers)
    set(POLYHEDRAL_TEMPLATE_MATCHING_BOOST_TARGET Boost::headers)
elseif(TARGET Boost::boost)
    set(POLYHEDRAL_TEMPLATE_MATCHING_BOOST_TARGET Boost::boost)
elseif(TARGET boost::headers)
    set(POLYHEDRAL_TEMPLATE_MATCHING_BOOST_TARGET boost::headers)
elseif(TARGET boost::boost)
    set(POLYHEDRAL_TEMPLATE_MATCHING_BOOST_TARGET boost::boost)
else()
    message(FATAL_ERROR "PolyhedralTemplateMatching requires a Boost headers target")
endif()

set(LIB_SOURCES
    ${CMAKE_SOURCE_DIR}/src/polyhedral_template_matching_service.cpp
    ${CMAKE_SOURCE_DIR}/src/polyhedral_template_matching.cpp
    ${CMAKE_SOURCE_DIR}/src/ptm_neighbor_finder.cpp
    ${CMAKE_SOURCE_DIR}/src/ptm_cluster_builder.cpp
    ${CMAKE_SOURCE_DIR}/src/structure_analysis_ptm.cpp
)

add_library(${PROJECT_NAME}_lib STATIC ${LIB_SOURCES})
set_target_properties(${PROJECT_NAME}_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(${PROJECT_NAME}_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/../CoreToolkit/dependencies/ptm
)
target_precompile_headers(${PROJECT_NAME}_lib PRIVATE <volt/core/volt.h>)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC
    ${POLYHEDRAL_TEMPLATE_MATCHING_BOOST_TARGET}
    TBB::tbb
    coretoolkit::coretoolkit
    structure-identification::structure-identification
    spdlog::spdlog
    nlohmann_json::nlohmann_json
)

add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_precompile_headers(${PROJECT_NAME} PRIVATE <volt/core/volt.h>)
target_link_libraries(${PROJECT_NAME} PRIVATE
    TBB::tbb
    ${PROJECT_NAME}_lib
)

if(POLYHEDRAL_TEMPLATE_MATCHING_REFRESH_LOCAL_DEPS)
    set(_ptm_local_dep_commands)

    if(EXISTS "${POLYHEDRAL_TEMPLATE_MATCHING_LOCAL_CORETOOLKIT_BUILD_DIR}/CMakeCache.txt")
        list(APPEND _ptm_local_dep_commands
            COMMAND "${CMAKE_COMMAND}" --build "${POLYHEDRAL_TEMPLATE_MATCHING_LOCAL_CORETOOLKIT_BUILD_DIR}" --parallel
        )
    endif()

    if(EXISTS "${POLYHEDRAL_TEMPLATE_MATCHING_LOCAL_STRUCTURE_IDENTIFICATION_BUILD_DIR}/CMakeCache.txt")
        list(APPEND _ptm_local_dep_commands
            COMMAND "${CMAKE_COMMAND}" --build "${POLYHEDRAL_TEMPLATE_MATCHING_LOCAL_STRUCTURE_IDENTIFICATION_BUILD_DIR}" --parallel
        )
    endif()

    if(_ptm_local_dep_commands)
        add_custom_target(${PROJECT_NAME}_local_deps
            ${_ptm_local_dep_commands}
            COMMENT "Refreshing local development dependencies for ${PROJECT_NAME}"
            VERBATIM
        )
        add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_local_deps)
    endif()
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(TARGETS ${PROJECT_NAME}_lib DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)
