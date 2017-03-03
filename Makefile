.DEFAULT: all

CC        = gcc
CFLAGS    = -O3 -Wno-unused-result -MMD
VENDOR := $(sort $(dir $(wildcard vendor/*/)))
SRC_C := $(wildcard src/*.c)
SRC_O := $(SRC_C:.c=.o)

EXE_C := $(wildcard samples/*.c)
EXE   := $(EXE_C:.c=)
INCLUDE   = -Isrc/ $(addprefix -I,$(VENDOR))

src/%.h: src/%.c
	@touch $@

src/%.o: src/%.c src/%.h
	@echo "Compiling $@"
	@$(CC) -o $@ $(CFLAGS) -c $<

samples/%: samples/%.c
	@echo "Compiling $@"
	@$(CC) -o $@ $(CFLAGS) $(INCLUDE) $< $(shell $(CC) -MM $(INCLUDE) $<|sed 's/\\//;s/\s/\n/g;s/\.h/.o/g'|tail -n +3|sort|uniq|xargs) -lm 

all: $(SRC_O) $(EXE)

clean:
	rm -f $(SRC_O) $(EXE) $(SRC_O:.o=.d) $(EXE_C:.c=.d)

vendor-all:
	cd vendor/hungarian && make all

vendor-clean:
	cd vendor/hungarian && make clean

# -include $(SRC_O:.o=.d)
-include $(EXE_C:.c=.d)
# DO NOT DELETE
