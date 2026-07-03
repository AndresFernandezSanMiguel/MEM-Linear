cmake_minimum_required(VERSION 3.14)
project(MEM_saoas LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# ---------------------------------------------------------------------------
# Dependencias: Eigen y Spectra (ambas header-only), descargadas
# automaticamente en tiempo de configuracion via FetchContent.
# ---------------------------------------------------------------------------
include(FetchContent)

FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
)
# No necesitamos compilar los tests/demos de Eigen
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
  spectra
  GIT_REPOSITORY https://github.com/yixuan/spectra.git
  GIT_TAG v1.0.1
)

FetchContent_MakeAvailable(eigen spectra)

# ---------------------------------------------------------------------------
# Ejecutable principal
# ---------------------------------------------------------------------------
add_executable(saoas src/saoas.cpp)

target_include_directories(saoas PRIVATE
  ${spectra_SOURCE_DIR}/include
)

target_link_libraries(saoas PRIVATE Eigen3::Eigen)

target_compile_options(saoas PRIVATE -O3 -Wall)
