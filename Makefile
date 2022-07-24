PROGRAM=ss
CC=gcc
LFLAGS=-lflint -lmpfr -lgmp -lpthread
CFLAGS=-Wall -O3

run: $(PROGRAM)
		./$(PROGRAM) 

$(PROGRAM): $(PROGRAM).o
		$(CC) $(PROGRAM).o -o $(PROGRAM) $(LFLAGS)

$(PROGRAM).o: $(PROGRAM).c
		$(CC) -c $(PROGRAM).c $(CFLAGS)

clean:
		rm *.o $(PROGRAM)
