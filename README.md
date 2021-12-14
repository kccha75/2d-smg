# Spectral Multigrid

Spectral Multigrid code written in Matlab.

## About:
Solves one-dimensional or two-dimensional PDEs using Spectral Multigrid.

Each direction can be discretised by either Chebyshev or Fourier collocation methods.

Supports Dirichlet, Neumann and mixed boundary conditions for Chebyshev.

## Important Notes:

Code is separated into each multigrid components:

> **Multigrid components**
* **Relaxation** - Minimum residual relaxation
* **Prolongation** - Fourier and Chebyshev
* **Restriction** - Fourier and Chebyshev
* **Preconditioner** - 2nd-order two-dimensional finite difference operator
* **Solver** - Conjugate gradient / Bicgstab, two-dimensional exact spectral matrix inversion
* **Scheme** - two-dimensional Poisson's equation

Prolongation and Restriction can easily be generalised for higher dimensions or finite difference scheme

Preconditioner, Solver and Scheme can easily be extended for different PDEs and higher dimensions or finite difference scheme

Multigrid, Minimum residual relaxation and conjugate gradient / bicgstab apply to any problems and dimensions provided the correct scheme is inputted. 

## Example solution
Poisson's equation with y periodic boundary conditions, x=-1 homogeneous Dirichlet, x=1 homogeneous Neumann

Exact solution is given by $$(cosh(1/2*(x-1))-cosh(1)).*exp(sin(x))$$

![Plot of numerical solution](./Images/plot1.pdf)

![Contours of numerical solution](./Images/plot2.pdf)

![Numerical errors](./Images/plot3.pdf)
