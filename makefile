all:
	g++ -lgsl -o relaxor -g main.cpp sistema.cpp impresor.cpp


relaxor:	main.o sistema.o impresor.o
	g++ main.o sistema.o impresor.o -o relaxor

main.o: main.cpp
	g++ -c main.cpp

sistema.o: sistema.cpp
	g++ sistema.cpp

impresor.o: impresor.cpp
	g++ impresor.cpp

clean:
	rm -rf *o relaxor
