default: debug

debug: debug-bint.o bigfact.c
	gcc -g -Wall bigfact.c debug-bint.o -o bigfact -lpthread

release: bint.o bigfact.c
	gcc -O3 -fno-inline-functions bigfact.c bint.o -o bigfact -lpthread

%-test: test/%.c debug-bint.o
	gcc -g -Wall $< debug-bint.o -o $@ -lpthread

debug-%.o: %.c
	gcc -g -Wall -c $< -o $@ -lpthread

%.o: %.c
	gcc -O3 -fno-inline-functions -c $< -o $@ -lpthread

clean:
	rm *.o
