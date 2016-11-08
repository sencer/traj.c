.DEFAULT: all
export CC        = gcc
export CFLAGS    =

all:   CFLAGS += -O3 -Wno-unused-result
debug: CFLAGS += -g -Wall -pedantic -fno-omit-frame-pointer -fsanitize=address

.PHONY: all debug clean

EXEC      = bond overlap unwrap nonbonded findbox travel wrap mass test h1h2 protonated vacancy sphere_vacancy
DEPS      = bonding.o boxes.o crystal.o util.o lammpstrj.o xyz.o fragments.o hungarian/hungarian.o periodic_table.o

all: src-all vendor-all samples-all 
debug: src-debug samples-debug
clean: src-clean vendor-clean samples-clean 

src-all:
	cd src && make all

src-debug:
	cd src && make debug

src-clean:
	cd src && make clean

samples-all: src-all
	cd samples && make all

samples-clean:
	cd samples && make clean

vendor-all:
	cd vendor/hungarian && make all

vendor-clean:
	cd vendor/hungarian && make clean
