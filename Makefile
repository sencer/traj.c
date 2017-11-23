.DEFAULT: all

CC        = gcc
CFLAGS    = -O3 -Wno-unused-result -MMD -fsanitize=address
VENDOR := $(sort $(dir $(wildcard vendor/*/)))
SRC_C := $(wildcard src/*.c)
SRC_O := $(SRC_C:.c=.o)

EXE_C := $(wildcard samples/*.c)
EXE   := $(EXE_C:.c=)
INCLUDE   = -Isrc/ $(addprefix -I,$(VENDOR))

all: vendor-all $(SRC_O) $(EXE)

src/%.h: src/%.c
	@touch $@

src/%.o: src/%.c src/%.h
	@echo "Compiling $@"
	@$(CC) -o $@ $(CFLAGS) -c $<

samples/%: samples/%.c
	@echo "Compiling $@"
	@$(CC) -o $@ $(CFLAGS) -fopenmp $(INCLUDE) $< $(shell $(CC) -MM $(INCLUDE) $<|sed 's/\\//;s/\s/\n/g;s/\.h/.o/g'|tail -n +3|sort|uniq|xargs) -lm 

samples/obrings: samples/obrings.cpp
	g++ -o samples/obrings -O3 -Isrc/ -I/usr/local/include/openbabel-2.0 $< -Lsrc/ -ltraj -lopenbabel -lm

clean: vendor-clean
	rm -f $(SRC_O) $(EXE) $(SRC_O:.o=.d) $(EXE_C:.c=.d)

vendor-all:
	cd vendor/hungarian && make all

vendor-clean:
	cd vendor/hungarian && make clean

# -include $(SRC_O:.o=.d)
-include $(EXE_C:.c=.d)
# DO NOT DELETE
