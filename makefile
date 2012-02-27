#Opciones variables
CFLAGS=-Wall -lgsl -c


all:	relaxor
#	g++ -lgsl -o relaxor -g main.cpp sistema.cpp impresor.cpp


relaxor:	main.o sistema.o impresor.o
	g++ main.o sistema.o impresor.o -o relaxor

main.o: main.cpp
	g++ $(CFLAGS) main.cpp

sistema.o: sistema.cpp
	g++ $(CFLAGS) sistema.cpp

impresor.o: impresor.cpp
	g++ $(CFLAGS) impresor.cpp

clean:
	rm -rf *o relaxor
