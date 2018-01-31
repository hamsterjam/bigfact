default: debug

debug:
	gcc -g -Wall -c bint.c -o bint.o -lpthread
	gcc -g -Wall bigfact.c bint.o -o bigfact -lpthread

release:
	gcc -O3 -fno-inline-functions -c bint.c -o bint.o -lpthread
	gcc -O3 -fno-inline-functions bigfact.c bint.o -o bigfact -lpthread

clean:
	rm bint.o
