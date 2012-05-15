#Opciones variables
EXE=relaxor
SOURCES=sistema.cpp impresor.cpp
OBJECTS=$(SOURCES:.cpp=.o)
CFLAGS= -march=native -O2 -pipe -Wall -g -lgsl -lgslcblas

all: $(EXE)

$(EXE): main.o $(OBJECTS)
	g++ $^ -o $@ $(CFLAGS)

tests: tester.o $(OBJECTS)
	g++ $^ -o $@ $(CFLAGS)

.cpp.o:
	g++ -c $< $(CFLAGS)

clean:
	rm -rf *o *~ $(EXE) *pyc


