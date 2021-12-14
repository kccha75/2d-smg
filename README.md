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
Poisson's equation in the rectangular domain

<img src="https://latex.codecogs.com/svg.image?x\in&space;[-1,1],\quad&space;y\in[-\pi,\pi]" title="x\in [-1,1],\quad y\in[-\pi,\pi]" />

* Periodic boundary conditions in <img src="https://latex.codecogs.com/svg.image?y" title="y" />
* Homogeneous Dirichlet on <img src="https://latex.codecogs.com/svg.image?x=-1" title="x=-1" />
* Homogeneous Neumann on <img src="https://latex.codecogs.com/svg.image?x=1" title="x=1" /> 

Exact solution is given by 

<img src="https://latex.codecogs.com/svg.image?u(x,y)=\left(\hbox{cosh}\left(\frac{1}{2}(x-1)\right)-\hbox{cosh}(1)\right)e^{\sin(y)}" title="u(x,y)=\left(\hbox{cosh}\left(\frac{1}{2}(x-1)\right)-\hbox{cosh}(1)\right)e^{\sin(y)}" />

## Numerical solution
![Plot of numerical solution](./Images/plot1.png)

## Contours of numerical solution
![Contours of numerical solution](./Images/plot2.png)

## Errors (compared to exact solution)
![Numerical errors](./Images/plot3.png)
