PROJECT(relaxor)
SET(CMAKE_CXX_FLAGS "-Wall -lgsl")
ADD_EXECUTABLE(relaxor impresor.cpp sistema.cpp main.cpp)
