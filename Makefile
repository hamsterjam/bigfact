default: debug

debug:
	gcc -g bigfact.c -o bigfact -lpthread

release:
	gcc -O3 bigfact.c -o bigfact -lpthread
