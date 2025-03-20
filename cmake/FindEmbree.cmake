find_path(EMBREE_INCLUDE_PATH embree4/rtcore.h
  ${CMAKE_SOURCE_DIR}/embree/include
  /usr/include
  /usr/local/include
  /opt/local/include)

if (APPLE)
find_library(EMBREE_LIBRARY NAMES embree4 PATHS
  ${CMAKE_SOURCE_DIR}/embree/lib-macos
  /usr/lib
  /usr/local/lib
  /opt/local/lib)
elseif (WIN32)
find_library(EMBREE_LIBRARY NAMES embree4 PATHS
  ${CMAKE_SOURCE_DIR}/embree/lib-win32)
else ()
find_library(EMBREE_LIBRARY NAMES embree4 PATHS
  ${CMAKE_SOURCE_DIR}/embree/lib-linux
  /usr/lib
  /usr/local/lib
  /opt/local/lib)
endif ()

if (EMBREE_INCLUDE_PATH AND EMBREE_LIBRARY)
  set(EMBREE_FOUND TRUE)
endif ()


# Locate FFTW include directory
find_path(FFTW_INCLUDE_DIR fftw3.h HINTS /opt/homebrew/include /usr/local/include)

# Locate FFTW library
find_library(FFTW_LIBRARY NAMES fftw3 HINTS /opt/homebrew/lib /usr/local/lib)

# Check if FFTW was found
if (FFTW_INCLUDE_DIR AND FFTW_LIBRARY)
    message(STATUS "Found FFTW: ${FFTW_LIBRARY}")
else()
    message(FATAL_ERROR "FFTW not found. Set FFTW_INCLUDE_DIR and FFTW_LIBRARY manually.")
endif()

# Add the include directory
include_directories(${FFTW_INCLUDE_DIR})


