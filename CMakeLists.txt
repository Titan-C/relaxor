cmake_minimum_required(VERSION 2.8)

project(relaxor)
execute_process(COMMAND pkg-config --libs gsl OUTPUT_VARIABLE LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

FIND_PACKAGE(Boost 1.46.1)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}" "/usr/include/python2.7")
  SET(Boost_USE_STATIC_LIBS OFF)
  SET(Boost_USE_MULTITHREADED ON)
  SET(Boost_USE_STATIC_RUNTIME OFF)
  FIND_PACKAGE(Boost 1.45.0 COMPONENTS python)

  ADD_LIBRARY(mat SHARED libmat.cpp)
  TARGET_LINK_LIBRARIES(mat ${Boost_LIBRARIES} ${LIBS})
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
ENDIF()