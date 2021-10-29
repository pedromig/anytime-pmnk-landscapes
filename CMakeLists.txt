cmake_minimum_required(VERSION 3.11)

if(${CMAKE_VERSION} VERSION_LESS 3.11)
  cmake_policy(VERSION ${CMAKE_MINOR_VERSION}.${CMAKE_MAJOR_VERSION})
endif()

# Define project properties
project(
  APMNKL
  VERSION 0.2.0
  LANGUAGES CXX
  DESCRIPTION "Library providing implementations of search \
               heuristics for the pmnk-landscapes problem."
)

# Turn off C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# Create compile commands json for clang tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# FetchContent_MakeAvailable implementation for cmake versions prior to 3.14
include(FetchContent)
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  macro(FetchContent_MakeAvailable NAME)
    fetchcontent_getproperties(${NAME})
    if(NOT ${NAME}_POPULATED)
      fetchcontent_populate(${NAME})
      add_subdirectory(${${NAME}_SOURCE_DIR} ${${NAME}_BINARY_DIR})
    endif()
  endmacro()
endif()

# Get target dependency (CLI11)
fetchcontent_declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.1.2
)
fetchcontent_makeavailable(CLI11)

# Library documentation
option(APMNKL_BUILD_DOCS "Build documentation" ON)

if(CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("-- Setting up documentation build")

    add_custom_target(
      doc_doxygen ALL
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "-- Generating documentation target"
      VERBATIM
    )
  else(DOXYGEN_FOUND)
    message("-- Could not find doxygen executable")
  endif(DOXYGEN_FOUND)
endif()

# Compiler Warnings (https://github.com/lefticus/cppbestpractices/)
option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" ON)

set(
  MSVC_WARNINGS
  /W4
  /w14242
  /w14254
  /w14263
  /w14265
  /w14287
  /we4289
  /w14296
  /w14311
  /w14545
  /w14546
  /w14547
  /w14549
  /w14555
  /w14619
  /w14640
  /w14826
  /w14905
  /w14906
  /w14928
  /permissive-
)

set(
  CLANG_WARNINGS
  -Wall
  -Wextra
  -Wshadow
  -Wnon-virtual-dtor
  -Wold-style-cast
  -Wcast-align
  -Wunused
  -Woverloaded-virtual
  -Wpedantic
  -Wconversion
  -Wsign-conversion
  -Wnull-dereference
  -Wdouble-promotion
  -Wformat=2
  -Wimplicit-fallthrough
)

if(WARNINGS_AS_ERRORS)
  set(CLANG_WARNINGS ${CLANG_WARNINGS} -Werror)
  set(MSVC_WARNINGS ${MSVC_WARNINGS} /WX)
endif()

set(
  GCC_WARNINGS
  ${CLANG_WARNINGS}
  -Wmisleading-indentation
  -Wduplicated-cond
  -Wduplicated-branches
  -Wlogical-op
  #-Wuseless-cast
)

if(MSVC)
  set(PROJECT_WARNINGS ${MSVC_WARNINGS})
elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  set(PROJECT_WARNINGS ${CLANG_WARNINGS})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(PROJECT_WARNINGS ${GCC_WARNINGS})
else()
  message("No compiler warnings set for '${CMAKE_CXX_COMPILER_ID}' compiler.")
endif()

# Anytime pmnk-landscapes library
set(APMNKL-LIB apmnkl)

add_library(${APMNKL-LIB} INTERFACE)
target_compile_features(${APMNKL-LIB} INTERFACE cxx_std_17)
target_compile_options(${APMNKL-LIB} INTERFACE ${PROJECT_WARNINGS})

target_include_directories(
  ${APMNKL-LIB} INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Install directories for library targets
include(GNUInstallDirs)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/include" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# The library is an install target
install(TARGETS ${APMNKL-LIB} EXPORT ${APMNKL-LIB}-targets)

install(
  EXPORT ${APMNKL-LIB}-targets
  FILE ${APMNKL-LIB}-config.cmake
  NAMESPACE ${APMNKL-LIB}::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${APMNKL-LIB}"
)

# Use find_package on the installed package
export(
  TARGETS ${APMNKL-LIB}
  NAMESPACE ${APMNKL-LIB}::
  FILE ${APMNKL-LIB}-targets.cmake
)

# Register in the user cmake package registry
export(PACKAGE ${APMNKL-LIB})

# Anytime pmnk-landscapes sample app (app.cpp) for testing implemented features.
set(APP anytime-pmnk-landscapes)

add_executable(${APP} "${PROJECT_SOURCE_DIR}/apps/app.cpp")
target_compile_features(${APP} PRIVATE cxx_std_17)
target_compile_options(${APP} PRIVATE ${PROJECT_WARNINGS})

target_link_libraries(${APP} PRIVATE ${APMNKL-LIB})
target_link_libraries(${APP} PRIVATE CLI11::CLI11)
