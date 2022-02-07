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

all:	test_interface2d.x test_raytrace.x test_tomography.x

test_interface2d.x:	test_interface2d.c ../../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf
	./$@ && rm $@

test_raytrace.x:	test_raytrace.c ../../raytrace.c ../../grid2.c ../../atela.c ../../dynamic.c ../../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm
	./$@ && rm $@

test_tomography.x:	test_tomography.c ../../tomography.c ../../raytrace.c ../../grid2.c ../../atela.c ../../dynamic.c ../../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm
	./$@ && rm $@

clean:
	rm *.x

help:
	@echo "Use 'make' to run the tests"
