#Opciones variables
TEST=test
SOURCES=material.cpp
HEADERS=$(SOURCES:.cpp=.h)
OBJECTS=$(SOURCES:.cpp=.o)
LIBS=$(shell pkg-config --libs gsl)
CFLAGS= -g -march=native -O2 -pipe -Wall

$(TEST): tester.o $(OBJECTS)
	g++ $^ -o $@ $(CFLAGS) $(LIBS)
	./$@

.cpp.o: $(HEADERS)
	g++ -c $< $(CFLAGS)

clean:
	rm -vf *.o *~ $(TEST) *pyc *dat

run: $(EXE)
	python main.py
	python plotter.py

