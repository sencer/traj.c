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
	@$(CC) -o $@ $(DEPFLAGS) $(INCLUDE) $(CFLAGS) $@.c $(shell sed -n 'H;/^$$/{g;s/\\\?\s\+/\n/g;s/\.h/.o/g;s/$@:\|$@.c//g;p;q}' $(DEPDIR)/$@.d|sort|uniq|xargs) -lm
	@$(POSTCOMPILE)

clean:
	rm -f $(EXEC)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRC)))
