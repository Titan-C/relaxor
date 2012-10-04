#Opciones variables
EXE=relaxor
TEST=test
SOURCES=material.cpp experiment.cpp impresor.cpp
OBJECTS=$(SOURCES:.cpp=.o)
LIBS=$(shell pkg-config --libs gsl)
CFLAGS= -g -march=native -O2 -pipe -Wall

all: $(EXE)

$(EXE): main.o $(OBJECTS)
	g++ $^ -o $@ $(CFLAGS) $(LIBS)

$(TEST): tester.o $(OBJECTS)
	g++ $^ -o $@ $(CFLAGS) $(LIBS)
	./$@

.cpp.o:
	g++ -c $< $(CFLAGS)

clean:
	rm -vf *.o *.*~ $(EXE) $(TEST) *pyc *dat

run: $(EXE)
	./$(EXE) 12 1 0.1 "9 0.4 1" 0.3 1000
	python plotter.py
	rm *dat
