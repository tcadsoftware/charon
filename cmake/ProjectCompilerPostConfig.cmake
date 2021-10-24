SET(Trilinos_ENABLE_Kokkos ${${PROJECT_NAME}_ENABLE_Kokkos} CACHE BOOL "Setting Trilinos_ENABLE_Kokkos to ${PROJECT_NAME}_ENABLE_Kokkos in ProjectCompilerPostConfig.cmake.")
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/../Trilinos/cmake/ProjectCompilerPostConfig.cmake)
