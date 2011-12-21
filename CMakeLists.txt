PROJECT(relaxor)
SET(CMAKE_CXX_FLAGS "-Wall -lgsl -march=native -O2 -pipe")
ADD_EXECUTABLE(relaxor impresor.cpp sistema.cpp main.cpp)
