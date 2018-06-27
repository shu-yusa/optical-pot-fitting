F90 = gfortran
F90omp = gfortran -xW -vec-report1 -openmp -openmp_report2
objs = $(.o, .mod)

RM := rm

a.out : global_constant.o  Coulomb_mod.o fit_opt_pot.o
	${F90omp} global_constant.o Coulomb_mod.o fit_opt_pot.o
global_constant.o : global_constant.f90
	${F90} -c global_constant.f90
Coulomb_mod.o : Coulomb.f90
	${F90omp} -c Coulomb_mod.f90
fit_opt_pot.o : global_constant.f90 Coulomb_mod.f90 fit_opt_pot.f90
	${F90omp} -c fit_opt_pot.f90

.PHONY : clean
clean : 
	$(RM) *.o *.mod a.out gnu* *.eps
