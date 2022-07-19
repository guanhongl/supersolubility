PROGRAM=ss
PROGRAM_T=ss_threaded
CC=gcc
LFLAGS=-lflint -lmpfr -lgmp -lpthread
CFLAGS=-Wall -O3

all: run $(PROGRAM_T)

run: $(PROGRAM)
		./$(PROGRAM) 

$(PROGRAM): $(PROGRAM).o
		$(CC) $(PROGRAM).o -o $(PROGRAM) $(LFLAGS)

$(PROGRAM).o: $(PROGRAM).c
		$(CC) -c $(PROGRAM).c $(CFLAGS)

$(PROGRAM_T): $(PROGRAM_T).o
		$(CC) $(PROGRAM_T).o -o $(PROGRAM_T) $(LFLAGS)

$(PROGRAM_T).o: $(PROGRAM_T).c
		$(CC) -c $(PROGRAM_T).c $(CFLAGS)

clean:
		rm *.o $(PROGRAM) $(PROGRAM_T)
