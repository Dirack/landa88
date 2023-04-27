# Makefile
# 
# Purpose: Compile and run unit tests.
# 
# Site: https://dirack.github.io
# 
# Versions 1.0
# 
# Programmer: Rodolfo A C Neves 25/12/2021
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

SRC := tomography.c raytrace.c grid2.c  atela.c dynamic.c interface2d.c layer2d.c model2d.c
TESTS := $(patsubst %.c,%.x,$(wildcard test_*.c))
LANDA88 := $(RSFSRC)/user/landa88
LIB := atestlib
UNITY := Unity/unity.c
GCC := gcc
GCC_OPT := $(UNITY) -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L. -L$$RSFSRC/lib -lrsf -lm -l$(LIB)
DEP := $(addprefix $(LANDA88)/, $(SRC))

# Function to compile and run a test
# The test executable will be removed after running
define RUN_TEST
	$(GCC) $^ $(GCC_OPT) -o $@
	./$@; rm $@
endef

all:	$(TESTS)

# Create a static lib with dependencies
lib$(LIB).a:	$(DEP)
	gcc -c $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -lm -lrsf
	ar -rc $@ $(patsubst %.c,%.o,$^)
	ranlib $@

# Rule for running a test
%.x:	%.c lib$(LIB).a
	$(call RUN_TEST)

clean:
	rm -f *.o *.x lib$(LIB).a

help:
	@echo "Use 'make' to run all the tests"
	@echo "Or run 'make test_filename.x' to run a specific test called filename"
	@echo "Tests available: $(TESTS)"

.PHONY: clean help
