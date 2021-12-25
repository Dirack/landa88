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

all:	test_interface2d.x
	./$< && rm $<

%.x:	%.c ../../interface2d.c ../../cubicsplineint.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf

clean:
	rm *.x

help:
	@echo "Use 'make' to run the tests"
