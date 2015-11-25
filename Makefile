CC=gcc
CFLAGS= -lm -Wall -pedantic -std=c99
EXEC = bond
DEPS = bonding.o boxes.o crystal.o util.o lammpstrj.o

.PHONY: all debug clean

all: CFLAGS += -fPIC -O3
all: $(EXEC)

debug: CFLAGS += -g
debug: $(EXEC)

$(EXEC): % : %.c libbond.a
	$(CC) $@.c -L. -lbond $(CFLAGS) -o $@

libbond.a:$(DEPS)
	ar rvs libbond.a $(DEPS)

%.o:%.c %.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(EXEC) libbond.a

%.c: ;
%.h: ;
Makefile: ;
