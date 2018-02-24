cmake_minimum_required(VERSION 2.8)

project(inti)

if (USE_TCL_PLUGIN)
    find_package(Tcl)
    if(TCL_FOUND)
        include_directories(${TCL_INCLUDE_PATH})
        set(EXTRA_LIBS ${EXTRA_LIBS} ${TCL_LIBRARY})
        set(ALL_FILES ${ALL_FILES} plugin-tcl.c plugin-tcl.h)
        set(_TCL_VERSION_ ${TCL_VERSION})
        set(_TCL_PLUGIN 1)
    endif(TCL_FOUND)
endif(USE_TCL_PLUGIN)

find_package(TCL REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/include
        ${TK_INCLUDE_PATH}
        ${TCL_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "-DNDEBUG -O3 ${CMAKE_CXX_FLAGS}")

SET(CMAKE_VERBOSE_MAKEFILE ON)
add_definitions(-DUSE_NON_CONST -DUSE_COMPOSITELESS_PHOTO_PUT_BLOCK)
add_executable(intiwish src/photo.cc src/tkMain.c src/morfo.cc src/contour.cc
        src/feature.cc src/filtre.cc src/histo.cc src/Image.cc
        src/hough_line.cc src/hough_circle.cc src/hough_track.cc)
add_custom_command(TARGET intiwish POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy intiwish ${CMAKE_SOURCE_DIR}/.)

target_link_libraries(intiwish  ${TCL_LIBRARY} ${TK_LIBRARY})