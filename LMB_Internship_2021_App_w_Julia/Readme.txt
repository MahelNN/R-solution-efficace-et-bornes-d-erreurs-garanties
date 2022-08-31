######################################################################
## @Autor: Mahel Ndoumba N.                                         ##  
## LMB internship 2021 : Efficient resolution and guaranteed bounds ##
##                         for eigenvalue problems in electronic    ##
##                                structure calculations            ##
######################################################################

The files in _h.jl contain the subroutines used in _.jl files.

00Exact_Sol_h.jl
00Exact_Sol.jl
Compute the parameters for the analytical solution and some plots of the solution

01FEM_h.jl
01FEM.jl
Implement the finite element method, some plots for the numerical solution and some error analysis

02GalMth_GBases_h.jl
02GalMth_GBases.jl        (Optimization of \sigma with respect the approximated eigenvalue)
02GalMth_GBases_res.jl     (Optimization of \sigma with respect the dual norm of the residual)
Implement the Galerkin method with Gaussian bases and some plots for the numerical solution and some error analysis