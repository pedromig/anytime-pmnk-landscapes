cmake_minimum_required(VERSION 3.11...3.16.3)

if(NOT DEFINED PROJECT_NAME)
  set(CMAKE_CXX_EXTENSIONS OFF)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

project(
  anytime-pmnk-landscapes
  VERSION 1.0
  LANGUAGES CXX
  DESCRIPTION "Anytime analysis of algorithms for the pmnk-landscapes problem"
)

if(${CMAKE_VERSION} VERSION_LESS 3.11)
  cmake_policy(VERSION ${CMAKE_MINOR_VERSION}.${CMAKE_MAJOR_VERSION})
endif()

include(FetchContent)
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  macro(FetchContent_MakeAvailable NAME)
    fetchcontent_getproperties(${NAME})
    string(TOLOWER "${NAME}" name)
    if(NOT ${name}_POPULATED)
      fetchcontent_populate(${name})
      add_subdirectory(${${name}_SOURCE_DIR} ${${name}_BINARY_DIR})
    endif()
  endmacro()
endif()


option(BUILD_DOC "Build documentation" ON)

if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")

  find_package(Doxygen)
  if (DOXYGEN_FOUND)
      set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
      set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

      configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
      message("-- Documentation build started")

      add_custom_target(doc_doxygen ALL
          COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
          COMMENT "-- Generating documentation..."
          VERBATIM )
  else (DOXYGEN_FOUND)
    message("-- Could not find doxygen executable...")
  endif(DOXYGEN_FOUND)
endif()

file(
  GLOB_RECURSE PMNK CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/src/*.cpp"
  "${PROJECT_SOURCE_DIR}/src/*.hpp"
)

add_executable(${PROJECT_NAME} ${PMNK})

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic -pedantic-errors -Wconversion -Wsign-conversion)

# fetchcontent_declare(
#   Catch2
#   GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#   GIT_TAG v2.13.1
# )

# fetchcontent_makeavailable(Catch2)
# target_link_libraries(${PROJECT_NAME} PRIVATE Catch2::Catch2)

fetchcontent_declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v1.9.1
)

fetchcontent_makeavailable(CLI11)
target_link_libraries(${PROJECT_NAME} PRIVATE CLI11::CLI11)
