# Quasi-3D stationary solutions of an incompressible flow in a closed basin with wind forcing
## _PHISICAL FORMULATION_

In order to simulate the behavior of a fluid in a close basin we should start from the navier stokes equations:

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/1.gif)

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/2.gif)

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/3.gif)


Fx, Fy, Fz in general are the resultant of any external force acting on our system:
in our case we have the gravity force acting along the z-axis 

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/4.gif)

and the coriolis inertial force acting along x and y axes as result of earth rotation (the centrifugal force is neglectable)

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/5.gif)

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/6.gif)

f is the coriolis parameter, is about 10E-4 rad/s mutiplied by the sine of the latitude.
For a wide shallow basin we can simplify the third equation with the hidrostatic approximation:

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/7.gif)

Considering a uniform water density we can write the spatial derivatives of p, for example ∂p∂x as

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/8.gif)

with n the elevation of the surface of the basin

Then we consider incompressibility

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/9.gif)

Integrating on all the water column we obtain

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/10.gif)

with U,V integrals of the velocity over all the water column 

The system of differential equations now looks like this

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/11.gif)

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/12.gif)

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/8.gif)

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/10.gif)

As boundary condition for a closed basin seems reasonable a no-slip condition and a null pressure gradient

![u_boundary](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/13.gif)

![v_boundary](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/14.gif)

![∂p/∂xi(boundary)](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/15.gif)



the system is closed and, given an initial condition, we are ready to numerically solve it.

## _NUMERICAL FORMULATION_

### _Definition of the grid_

The grids used are the easiest possible: uniform 3-D  solidal with our x,y,z axes.
The only complication is the using of two grids, one for the velocity and one staggered grid for the pressure, with each cell of the second one centered in the intesections of the cells of the first one.

For general purpose and definition of our algorithms, a Staggered CFD Grid shall be used. It ensures the presence of very real non-zero pressure gradient across the nodes in any condition, even in the case of a checkered grid. The staggering also ensures realistic behaviors of the descretized momentum equations for spatially oscillating pressures. Also, the direction of the velocity vectors are exact.

So the velocity grid will have n-horizontal per m-vertical points and the pressure grid will have n+1 per m+1.
Therefore the (i,j,k) cell of the velocity grid will be centered between the (i,j,k),(i+1,j,k),(i,j+1,k),(i+1,j+1,k) cells of the pressure grid.

### _Discretization of the equations_

The Pressure equation can be easily discetized, calling i',j',k' the x,y,z indices in the pressure grid and i , j , k the x,y,z indices in the velocity grid 
we can calculate U(i,j) and V(i,j) as

![U_integral](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/U.gif)

![V_integral](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/V.gif)



Since U,V lie on a different grid than the p points, we shall calculate the discrete derivatives ΔU/Δx(i',j') , ΔV/Δy(i',j') by mediating on the direction ortogonal to the derivative itself.

![du_dx](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/dU_dx.gif)

![dv_dy](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/dV_dy.gif)



Finally we can write our nu equation discretized

![nu](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/nu.gif)

The discretization of the NS equations runs in the same staggering problem encountered in the discretization of the basin-level-equation but we should convert our nu derivative terms from the pressure grid to the velocity one so we can write the derivatives:

![d_eta_dx](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/d_eta_dx.gif)

![d_eta_dy](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/d_eta_dy.gif)

Furthermore we should find a way to discretize the advection term to let each cell communicate to all the adiacent and conservate locally the momentum.
Notice that if we lazy discretize the derivative as done before

![udxu1](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/udxu1.gif)

each cell communicate only with the following one and not with the previous one, so the most obvious solution is the centered finite difference

![udxu2](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/udxu2.gif)

so we can write the N.S. equations as

![NS_u](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/NS_u.gif)

![NS_v](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/NS_v.gif)

With the coriolis force 

![F_co](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/Fco.gif)

And the wind stress, given the 3m wind field

![F_wind](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/F_wind.gif)


In order to find the stationary solutions we will calculate the relative variations of velocity and pressure each time step, then we will sum the square values 
and we will look for a solution with a small enough value of total relative variation.

![TRV](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/TRV.gif)

![Trv_min](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/Figures/Trv_min.gif)

## _STRUCTURE OF THE CODE_



## _HOW TO RUN THE CODE_

To run the code you have to provide seven text files, in the folder [initial_conditions](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions):

*[u_top.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/u_top.txt) : n*m matrix with the values of the initial zonal velocity of the upper layer(m/s)

*[u_bot.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/u_bot.txt) : n*m matrix with the values of the initial zonal velocity of the lower layer(m/s)

*[v_top.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/v_top.txt) : n*m matrix with the values of the initial meridional velocity of the upper layer(m/s)

*[v_bot.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/v_bot.txt) : n*m matrix with the values of the initial meridional velocity of the lower layer(m/s)

*[u_wind.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/u_wind.txt) : n*m matrix with the values of the zonal wind velocity

*[v_wind.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/v_wind.txt) : n*m matrix with the values of the meridional wind velocity

*[const.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/const.txt) : contains all the phisical costants, the grid dimension and the number of cells in each direction in the velocity grid (to verify that data are correctly privided)

The user must provide phisically adequate values of variables:

*Current velocity must be below 1[m/s]

*Wind velocity must be below 30[m/s]

*Surface elevation must be below 2[m] and its mean value must be null

*The grid itself must be at least 3*3*2

To verify the correctnes of the provided data the user should compile [input_test.py](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/input_test.py).

Then and only then he will lauch the simulation compiling [cfd.py](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/cfd.py).























