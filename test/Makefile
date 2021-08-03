# Makefile
# 
# Purpose: Run all tests with 'make'.
# 
# Site: https://dirack.github.io
# 
# Version 1.0
# 
# Programer: Rodolfo A C Neves (Dirack) 05/10/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

DEP= ../raytrace.c ../grid2.c ../atela.c
SRCLANDA = ../Mlanda88tomo.c ../tomography.c $(DEP) ../vfsacrsnh_lib.c

# model files: Stacked section and parameters cube from CRE stack
PARCUBE=model/parametersCube.rsf
STKSEC=model/stackedSection.rsf

# Extract parameters from stacked section
NP=$(shell (<landa88tomo/model/stackedSection.asc sfget n1 n2; echo "n1*n2") | bc)
N1=$(shell <landa88tomo/model/stackedSection.asc sfget n1 parform=n)
D1=$(shell <landa88tomo/model/stackedSection.asc sfget d1 parform=n)
O1=$(shell <landa88tomo/model/stackedSection.asc sfget o1 parform=n)
N2=$(shell <landa88tomo/model/stackedSection.asc sfget n2 parform=n)
D2=$(shell <landa88tomo/model/stackedSection.asc sfget d2 parform=n)
O2=$(shell <landa88tomo/model/stackedSection.asc sfget o2 parform=n)
PARAM=np=$(NP) n1=$(N1) d1=$(D1) o1=$(O1) n2=$(N2) d2=$(D2) o2=$(O2)

# Function to call gdb scripts
# 1 is the number of the script
# 2 are the parameters
# 3 is the file to debug
define GDB
	gdb -q --command=gdb_scripts/gdb_script$(1) $(2) --args $(3)
endef

all:	landa88tomo.x
	@echo "Running integration tests, run 'make help' for more details"
	$(call GDB,1,-ex "r" -ex "q",landa88tomo.x)
	
setup1:	landa88tomo/pick.txt landa88tomo/parametersCube.asc
	@echo "Setup input files for test..."
	cd landa88tomo && scons -f SConscript1 $(PARAM) && scons model_0.rsf && scons sz.rsf && scons ns=$(shell wc -l landa88tomo/pick.txt | cut -d" " -f1) ni=1 nipsources_1.rsf angles_1.rsf gz_1.rsf

# Version of the landa88tomo program for GDB (GNU debugger)
landa88tomo:	landa88tomo.x GDB1

landa88tomo.x:	setup1 $(SRCSTEREO)
	@echo "Program landa88tomo.x for GDB built"
	gcc -DGDB_DEBUG -g $(SRCLANDA) -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm

GDB1:
	@echo "Running GDB script for landa88tomo.x"
	$(call GDB,1,-ex "b 193" -ex "r",landa88tomo.x)

clean:
	rm *.x

.PHONY: clean landa88tomo.x

help:
	@echo "Use 'make' to run unit tests and intergration test"
	@echo "Use 'make progname' to generate a GDB version of progname for debug"
	@echo "Programs available:"
	@echo "\tlanda88tomo"

landa88tomo/pick.txt:
	$(MAKE) -C landa88tomo picking

landa88tomo/parametersCube.asc:
	$(MAKE) -C landa88tomo parametersCube.asc
