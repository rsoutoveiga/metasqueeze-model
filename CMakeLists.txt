cmake_minimum_required(VERSION 3.5)

# Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    #if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

project(metapop_cpp17 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
#target_compile_features(${TARGET_NAME} PUBLIC cxx_std_17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options("/W4" "$<$<CONFIG:RELEASE>:/O2>")
else()
    add_compile_options(
        "-Wall" 
        "-Wextra"  
        "-Wpedantic" 
        "-Wshadow" 
        "-Wnon-virtual-dtor" 
        "-Wold-style-cast" 
        "-Wcast-align" 
        "-Wunused" 
        "-Woverloaded-virtual" 
        "-Wsign-conversion" 
        "-Wnull-dereference" 
        "-Wdouble-promotion" 
        "-Wformat=2"
        "-Wswitch-enum"
        "-Wswitch"
        "$<$<CONFIG:DEBUG>:-g>"
        "$<$<CONFIG:RELEASE>:-O3>"
        #"-Wconversion" 
        #"-Werror"
        )
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        add_compile_options("-stdlib=libc++")
    endif()
    
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        add_compile_options(
            "-Wmisleading-indentation" 
            "-Wduplicated-cond" 
            "-Wlogical-op" 
            "-Wuseless-cast"
            )
    endif()
endif()

include_directories(${CMAKE_SOURCE_DIR})

file(GLOB SOURCES
  src/*.cpp
)

file(GLOB HEADERS
  src/*.h
)

## this place the executable file in the project folder (i.e. where is this cmakelists file)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(metapop ${SOURCES} ${HEADERS})
