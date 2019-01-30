cmake_minimum_required(VERSION 3.5)

project(${B_PROJECT_NAME} VERSION ${B_VERSION_MAJOR}.${B_VERSION_MINOR}.${B_VERSION_PATCH} LANGUAGES C)

set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-std=gnu11 -fPIC -Wall -Wpedantic -Wextra -Wno-unused-parameter -Wno-format-security -Wno-int-conversion -Wno-incompatible-pointer-types -pthread")
# -Wextra -Wno-unused-variable
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_C_FLAGS_MINSIZEREL "-Os -DNDEBUG")

if (${B_MACHINE} MATCHES "x86_64")
  if (${CMAKE_C_COMPILER} MATCHES "clang")
    set(CMAKE_C_FLAGS_RELEASE "-Ofast3 -msse4.2 -DNDEBUG")
  endif()
  if (${CMAKE_C_COMPILER} MATCHES "gcc")
    set(CMAKE_C_FLAGS_RELEASE "-Ofast -msse4.2 -DNDEBUG")
  endif()
endif()

if (${B_MACHINE} MATCHES "aarch64")
  if (${CMAKE_C_COMPILER} MATCHES "clang")
    set(CMAKE_C_FLAGS_RELEASE "-Ofast3 -DNDEBUG")
  endif()
  if (${CMAKE_C_COMPILER} MATCHES "gcc")
    set(CMAKE_C_FLAGS_RELEASE "-Ofast -march=armv8-a -mtune=thunderx -DNDEBUG") # -mtune=cortex-a53
  endif()
endif()

set(CMAKE_BUILD_TYPE Release) #Debug, Release, RelWithDebInfo, MinSizeRel
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE SOURCE_FILES ${PROJECT_SOURCE_DIR}/*.c)

include_directories(${PROJECT_INCLUDE_DIR})

if(DEFINED QHULL_PREFIX)
  include_directories(${QHULL_PREFIX}/include)
  link_directories(${QHULL_PREFIX}/lib)
endif()

add_library(${B_PROJECT_NAME}_core OBJECT ${SOURCE_FILES})

set(ZLOG_CONF_PATH "/etc/semel/zlog.conf" CACHE STRING "Path to zlog configuration file")
target_compile_definitions(${B_PROJECT_NAME}_core PRIVATE ZLOG_CONF_PATH="${ZLOG_CONF_PATH}")

add_library(${B_PROJECT_NAME}_s STATIC $<TARGET_OBJECTS:${B_PROJECT_NAME}_core>)
add_library(${B_PROJECT_NAME}_sh SHARED $<TARGET_OBJECTS:${B_PROJECT_NAME}_core>)

target_link_libraries(${B_PROJECT_NAME}_s zlog Judy flint gsl gslcblas gmp mpfr qhull_r pthread m)
target_link_libraries(${B_PROJECT_NAME}_sh zlog Judy flint gsl gslcblas gmp mpfr qhull_r pthread m)

set_target_properties(${B_PROJECT_NAME}_s PROPERTIES OUTPUT_NAME ${B_PROJECT_NAME})
set_target_properties(${B_PROJECT_NAME}_sh PROPERTIES OUTPUT_NAME ${B_PROJECT_NAME} SOVERSION ${PROJECT_VERSION})

install(TARGETS ${B_PROJECT_NAME}_s ARCHIVE DESTINATION lib)
install(TARGETS ${B_PROJECT_NAME}_sh LIBRARY DESTINATION lib)
install(FILES ${PROJECT_INCLUDE_DIR}/${B_PROJECT_NAME}.h DESTINATION include/${B_PROJECT_NAME})
