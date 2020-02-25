# --------------------------------------------------------------------------- #
#    Main CMake file for MCFClass                                             #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8...3.15)
# TODO: set proper version
project(MCFClass VERSION 1.0.0 LANGUAGES CXX)

# ----- Requirements -------------------------------------------------------- #
# This is to use the find modules we provide
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(CPLEX)

# ----- Library ------------------------------------------------------------- #
add_library(MCFClass)
target_compile_features(MCFClass PUBLIC cxx_std_17)

target_sources(MCFClass PRIVATE
               MCFSimplex/MCFSimplex.C
               RelaxIV/RelaxIV.C
               SPTree/SPTree.C)

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake), that must be different
# for local build and installation.
target_include_directories(
        MCFClass PUBLIC
        $<INSTALL_INTERFACE:include/MCFClass>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFClass>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFClone>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFSimplex>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OPTUtils>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/RelaxIV>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SPTree>)

if (CPLEX_FOUND)
    target_sources(MCFClass PRIVATE MCFCplex/MCFCplex.C)
    target_include_directories(
            MCFClass PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFCplex>)
    target_link_libraries(MCFClass PUBLIC CPLEX::Cplex)
endif ()

# Alias for testing
add_library(MCFClass::MCFClass ALIAS MCFClass)

# ----- Install instructions ------------------------------------------------ #
# For installing the library and all its CMake modules on the system.
# Not required for local builds (see below).

include(GNUInstallDirs)

# Install library
install(TARGETS MCFClass
        EXPORT MCFClass-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/MCFClass/MCFClass.h
        ${CMAKE_CURRENT_SOURCE_DIR}/MCFClone/MCFClone.h
        ${CMAKE_CURRENT_SOURCE_DIR}/MCFSimplex/MCFSimplex.h
        ${CMAKE_CURRENT_SOURCE_DIR}/OPTUtils/OPTUtils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/RelaxIV/RelaxIV.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SPTree/SPTree.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MCFClass)

if (CPLEX_FOUND)
    install(FILES
            ${CMAKE_CURRENT_SOURCE_DIR}/MCFCplex/MCFCplex.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MCFClass)
endif ()

# Install target import file
install(EXPORT MCFClass-targets
        FILE MCFClassTargets.cmake
        NAMESPACE MCFClass::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFClass)

# Generate package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate package configuration file
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MCFClassConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFClass)

# Install package version and configuration files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFClass)

# Shipping custom find modules should not be encouraged because third
# party developers should provide them.
# We ship FindCPLEX because CPLEX does not on its own.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFClass)

# ----- Export from the build tree ------------------------------------------ #
# This is needed so other projects can find the library in a local build.
# It does not affect installation.

if (CMAKE_EXPORT_PACKAGE_REGISTRY)

    export(EXPORT MCFClass-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/MCFClassTargets.cmake
           NAMESPACE MCFClass::)

    # Export also the find modules for the local config file to use it
    configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
                   ${CMAKE_CURRENT_BINARY_DIR}/FindCPLEX.cmake
                   COPYONLY)

    # Register package in user's package registry (~/.cmake/registry)
    export(PACKAGE MCFClass)
endif ()
