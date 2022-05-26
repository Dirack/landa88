# Makefile
# 
# Purpose: Run all tests with 'make'.
# 
# Site: https://dirack.github.io
# 
# Version 1.0
# 
# Programmer: Rodolfo A C Neves (Dirack) 05/10/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

DEP= ../raytrace.c ../grid2.c ../atela.c ../dynamic.c ../interface2d.c ../layer2d.c ../model2d.c
TESTS= test_setup.x test_forward.x test_datamis.x test_vfsa.x

all:	$(TESTS)

test_setup.x:	test_setup.c ../raytrace.c ../grid2.c ../atela.c ../dynamic.c ../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm
	./$@ && rm $@

test_forward.x:	test_forward.c ../grid2.c ../atela.c ../dynamic.c ../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm
	./$@ && rm $@

test_datamis.x:	test_datamis.c ../grid2.c ../atela.c ../dynamic.c ../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm
	./$@ && rm $@
	
test_vfsa.x:	test_vfsa.c ../grid2.c ../atela.c ../dynamic.c ../interface2d.c Unity/unity.c
	gcc $^ -I$$RSFSRC/include -I$$RSFSRC/user/landa88 -L$$RSFSRC/lib -o $@ -lrsf -lm
	./$@ && rm $@


.PHONY:	$(TESTS)
