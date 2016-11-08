SRC    = $(wildcard *.c)
OBJ    = $(SRC:.c=.o)

DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

all debug: $(OBJ)

%.o: %.c
%.o: %.c $(DEPDIR)/%.d
	@echo "Compiling $< to $@"
	@$(CC) $(DEPFLAGS) $(CFLAGS) -c $<
	@$(POSTCOMPILE)

clean:
	rm -f *.o

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRC)))
