#Opciones variables
CFLAGS=-Wall -c -g
EXE=relaxor

all: relaxor

relaxor:	main.o sistema.o impresor.o
	g++ -g -Wall -lgsl $^ -o $(EXE)

main.o: main.cpp
	g++ $(CFLAGS) main.cpp

sistema.o: sistema.cpp
	g++ $(CFLAGS) sistema.cpp

impresor.o: impresor.cpp
	g++ $(CFLAGS) impresor.cpp

clean:
	rm -rf *o *~ $(EXE)
