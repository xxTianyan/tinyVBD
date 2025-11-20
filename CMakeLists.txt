cmake_minimum_required(VERSION 3.31)
project(tinyVBD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(EXTERN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern")

message(STATUS "EXTERN_DIR = ${EXTERN_DIR}")

# -----------  Eigen --------------
set(EIGEN_DIR "${EXTERN_DIR}/eigen-3.4.0")
if(NOT EXISTS "${EIGEN_DIR}/Eigen/Core")
    message(FATAL_ERROR
            "Eigen not found in ${EIGEN_DIR}\n"
            "Expected layout:\n"
            "  ${EIGEN_DIR}/Eigen/\n"
            "  ${EIGEN_DIR}/unsupported/\n")
endif()
add_library(eigen_headers INTERFACE)
target_include_directories(eigen_headers INTERFACE "${EIGEN_DIR}")
add_library(Eigen3::Eigen ALIAS eigen_headers)

# -----------  Raylib --------------
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory("${EXTERN_DIR}/raylib-master")

# -----------  ImGui --------------
set(IMGUI_DIR "${EXTERN_DIR}/imgui-master")
foreach(_f IN ITEMS imgui.cpp imgui_demo.cpp imgui_draw.cpp imgui_tables.cpp imgui_widgets.cpp)
    if(NOT EXISTS "${IMGUI_DIR}/${_f}")
        message(FATAL_ERROR "Missing ${IMGUI_DIR}/${_f}. Check extern/imgui is a full clone of Dear ImGui.")
    endif()
endforeach()

add_library(imgui STATIC
        "${IMGUI_DIR}/imgui.cpp"
        "${IMGUI_DIR}/imgui_demo.cpp"
        "${IMGUI_DIR}/imgui_draw.cpp"
        "${IMGUI_DIR}/imgui_tables.cpp"
        "${IMGUI_DIR}/imgui_widgets.cpp"
)
target_include_directories(imgui PUBLIC "${IMGUI_DIR}")

# ---------- rlImGui（Raylib + ImGui 的桥接库） ----------
set(RLIMGUI_DIR "${EXTERN_DIR}/rlImGui-main")
if(NOT EXISTS "${RLIMGUI_DIR}/rlImGui.cpp")
    message(FATAL_ERROR "Missing ${RLIMGUI_DIR}/rlImGui.cpp. Put rlImGui sources under extern/rlImGui")
endif()

add_library(rlImGui STATIC
        "${RLIMGUI_DIR}/rlImGui.cpp"
)
target_include_directories(rlImGui PUBLIC "${RLIMGUI_DIR}" "${IMGUI_DIR}")
target_link_libraries(rlImGui PUBLIC raylib imgui)

# ---------- My executable ----------
add_executable(tinyVBD src/main.cpp
        src/Mesh.cpp
        src/World.cpp
        src/BroadPhase.cpp)

target_include_directories(tinyVBD
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(tinyVBD PRIVATE Eigen3::Eigen rlImGui imgui raylib)

if(APPLE)
    target_compile_definitions(tinyVBD PRIVATE PLATFORM_MACOS)
elseif(WIN32)
    target_compile_definitions(tinyVBD PRIVATE PLATFORM_WINDOWS)
elseif(UNIX)
    target_compile_definitions(tinyVBD PRIVATE PLATFORM_LINUX)
endif()

if(MSVC)
    target_compile_options(tinyVBD PRIVATE /W4 /permissive- /utf-8 /Zc:__cplusplus)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(tinyVBD PRIVATE -Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(tinyVBD PRIVATE -Wall -Wextra -Wpedantic)
endif()

