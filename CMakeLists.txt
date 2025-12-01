cmake_minimum_required(VERSION 3.10)
project(maps_remesh)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set(SOURCES
    main.cpp
    BaseMesh.cpp
)

set(HEADERS
    BaseMesh.h
    Maps.h
    vertex_remove_data.h
)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_link_libraries(${PROJECT_NAME}
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
