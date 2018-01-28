default: debug

debug:
	gcc -g bigfact.c -o bigfact -lpthread

release:
	gcc -O3 -fno-inline-functions bigfact.c -o bigfact -lpthread
