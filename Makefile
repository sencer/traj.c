# CC        = clang-3.6
CC        = gcc
CFLAGS    = -std=c99
EXEC      = bond overlap unwrap
DEPS      = bonding.o boxes.o crystal.o util.o lammpstrj.o xyz.o

.PHONY: all debug clean
all: CFLAGS += -O3 -Wno-unused-result
debug: CFLAGS += -g -Wall -pedantic -fno-omit-frame-pointer -fsanitize=address

all: $(EXEC)
debug: $(EXEC)

$(EXEC): % : %.c $(DEPS)
	$(CC) $(CFLAGS) $@.c $(DEPS) -lm -o $@

$.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(DEPS) $(EXEC)

%.c: ;
%.h: ;
Makefile: ;
