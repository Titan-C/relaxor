#Opciones variables
EXE=relaxor
TEST=test
SOURCES=sistema.cpp impresor.cpp
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
	rm -vf *o *~ $(EXE) $(TEST) *pyc


