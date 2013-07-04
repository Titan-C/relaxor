#Opciones variables
EXE=relaxor
TEST=test
SOURCES=material.cpp experiment.cpp impresor.cpp
HEADERS=$(SOURCES:.cpp=.h)
OBJECTS=$(SOURCES:.cpp=.o)
LIBS=$(shell pkg-config --libs gsl)
CFLAGS= -g -march=native -O2 -pipe -Wall

all: $(EXE)

$(EXE): main.o $(OBJECTS) $(HEADERS)
	g++ $^ -o $@ $(CFLAGS) $(LIBS)

$(TEST): tester.o $(OBJECTS)
	g++ $^ -o $@ $(CFLAGS) $(LIBS)
	./$@

.cpp.o: $(HEADERS)
	g++ -c $< $(CFLAGS)

clean:
	rm -vf *.o *~ $(EXE) $(TEST) *pyc *dat

run: $(EXE)
	./$(EXE) 12 4 0.2 "7 0.2 1" 0.3 200
	python plotter.py
#	rm *dat

syscheck: $(EXE)
	valgrind --track-origins=yes --leak-check=full ./$(EXE) 5 1 0.1 "9 0.4 1" 0.3 1000
	valgrind --tool=cachegrind ./$(EXE) 5 1 0.1 "9 0.4 1" 0.3 1000
	valgrind --tool=callgrind ./$(EXE) 5 1 0.1 "9 0.4 1" 0.3 1000

