SRC       = $(wildcard *.c)
EXEC      = $(SRC:.c=)
INCLUDE   = -I../src
# DEPS      =  $(wildcard ../src/*.o) ../vendor/hungarian/hungarian.o

DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

all debug: $(EXEC)

$(EXEC): %: %.c
$(EXEC): %: %.c $(DEPDIR)/%.d
	@echo "Compiling $@"
	@$(CC) -o $@ $(DEPFLAGS) $(INCLUDE) $(CFLAGS) $@.c $(shell $(CC) -MM -I../src $@.c|sed 's/\\//;s/\s/\n/g;s/\.h/.o/g'|tail -n +3|sort|uniq|xargs) -lm
	@$(POSTCOMPILE)

clean:
	rm -f $(EXEC)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRC)))
	# @$(CC) -o $@ $(DEPFLAGS) $(INCLUDE) $(CFLAGS) $@.c $(shell sed -n 'H;/^$$/{g;s/\\\?\s\+/\n/g;s/\.h/.o/g;s/$@:\|$@.c//g;p;q}' $(DEPDIR)/$@.init|sort|uniq|xargs) -lm
