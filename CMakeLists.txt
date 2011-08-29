PROJECT(relaxor)
SET(CMAKE_CXX_FLAGS "-Wall -lgsl -O2")
ADD_EXECUTABLE(relaxor impresor.cpp sistema.cpp main.cpp)
