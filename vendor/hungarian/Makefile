CC = gcc

CFLAGS = -O3 -Wall -I. -fPIC

all: hungarian.o

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f hungarian.o
