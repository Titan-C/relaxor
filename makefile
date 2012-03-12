#Opciones variables
EXE=relaxor
SOURCES=main.cpp sistema.cpp impresor.cpp
OBJECTS=$(SOURCES:.cpp=.o)
CFLAGS= -march=native -O2 -pipe -Wall -g -lgsl -lgslcblas

all: $(EXE)

$(EXE): $(OBJECTS)
	g++ $(CFLAGS) $^ -o $@

.cpp.o:
	g++ $(CFLAGS) -c $<

clean:
	rm -rf *o *~ $(EXE)
