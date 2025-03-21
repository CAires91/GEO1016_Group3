cmake_minimum_required(VERSION 3.12)

get_filename_component(MODULE_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(PROJECT_NAME "easy3d_${MODULE_NAME}")
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        point_cloud_io.h
        image_io.h
        surface_mesh_io.h
        resources.h
        )

set(${PROJECT_NAME}_SOURCES
        point_cloud_io.cpp
        image_io.cpp
        surface_mesh_io.cpp
        surface_mesh_io_obj.cpp
        resources.cpp
        )


add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS})

# The resources directory
target_compile_definitions(
        ${PROJECT_NAME}
        PRIVATE
        "EASY3D_RESOURCES_DIR=\"${EASY3D_ROOT}/resources\""
)


set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "easy3d")

target_include_directories(${PROJECT_NAME} PRIVATE ${EASY3D_INCLUDE_DIR})

set(LASTOOLS_INCLUDE_DIR ${EASY3D_THIRD_PARTY}/LAStools/LASzip/src ${EASY3D_THIRD_PARTY}/LAStools/LASlib/inc)
target_include_directories(${PROJECT_NAME} PRIVATE ${LASTOOLS_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} easy3d_core easy3d_util)

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _CRT_SECURE_NO_DEPRECATE)
endif ()
