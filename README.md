Poissonsolver

Integral equation solver for Poisson's equation in the plane.

Todo: Setup & Build, Running code




#gfortran -fPIC -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -std=legacy -c -w solve_fmm.f -o solve_fmm.o


# gfortran -w -o int2.so -shared -fPIC lbfmm2d.o chebrouts.o tables8.o tree_routs4.o tree_vol_coeffs_cheb.o ./common/prini_new.o ./common/legeexps.o ./common/chebexps.o ./common/legetens.o ./common/chebtens.o ./common/voltab2d.o solve_fmm.o -lblas -llapack

