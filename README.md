# Closed Basin quasi-3D simulations of incompressible unsteady flow in a closed basin with idrostatic approx

PHISICAL FORMULATION

In order to simulate the behavior of a fluid in a close basin we should start from the navier stokes equations:

1 ∂u∂t+u∂u∂x+v∂u∂y=−1ρ∂p∂x+ν(∂2u∂x2+∂2u∂y2)+Fx
2 ∂v∂t+u∂v∂x+v∂v∂y=−1ρ∂p∂y+ν(∂2v∂x2+∂2v∂y2)+Fy
3 ∂w∂t+u∂v∂x+v∂v∂y=−1ρ∂p∂y+ν(∂2v∂x2+∂2v∂y2)+Fz

Fx, Fy, Fz in general are the resultant of any external force acting on our system:
in our case we have the gravity force acting along the z-axis 
4 Fz=-ρg

and the coriolis inertial force acting along x and y axes as result of earth rotation (the centrifugal force is neglectable)

5 Fx=fv
6 Fy=-fu

f is the coriolis parameter, is about 10E-4 rad/s mutiplied by the sine of the latitude.
For a wide shallow basin we can simplify the third equation with the hidrostatic approximation:
7 p=-ρgh

Considerin a uniform water density we can write the spatial derivatives of p, for example ∂p∂x as
8 ∂p∂x=-ρg∂n∂x

with n the elevation of the surface of the basin

Then we consider incompressibility

9 ∂u∂x+∂v∂y+∂w∂z=0

Integrating on all the water column we obtain
10 -(∂U∂x+∂V∂y)=∂n∂t

with U,V integrals of the velocity over all the water column 

The system of differential equations now looks like this

11 ∂u∂t+u∂u∂x+v∂u∂y=−1ρ∂p∂x+ν(∂2u∂x2+∂2u∂y2)+fv
12 ∂v∂t+u∂v∂x+v∂v∂y=−1ρ∂p∂y+ν(∂2v∂x2+∂2v∂y2)-fu
8 ∂p∂xi=-ρg∂n∂xi
10 -(∂U∂x+∂V∂y)=∂n∂t

As boundary condition for a closed basin seems reasonable a no-slip condition and a null pressure gradient
13 u(boundary)=0
14 v(boundary)=0
15 ∂p∂xi(boundary)=0


the system is closed and ,given an initial condition, we are ready to numerically solve it.

NUMERICAL FORMULATION

DEFINITION OF THE GRID

The grids used are the easiest possible: uniform 3-D  solidal with our x,y,z axes.
The only complication is the using of two grids, one for the velocity and one staggered grid for the pressure, with each cell of the second one centered in the intesections of the cells of the first one.

For general purpose and definition of our algorithms, a Staggered CFD Grid shall be used. It ensures the presence of very real non-zero pressure gradient across the nodes in any condition, even in the case of a checkered grid. The staggering also ensures realistic behaviors of the descretized momentum equations for spatially oscillating pressures. Also, the direction of the velocity vectors are exact.

So the velocity grid will have n-horizontal per m-vertical points and the pressure grid will have n+1 per m+1.

DISCRETIZATION OF THE EQUATIONS

The Pressure equation can be easily discetized, calling i',j',k' the x,y,z indices in the pressure grid and i , j , k the x,y,z indices in the velocity grid 
we can calculate U(i,j) and V(i,j) as

U
V

Since U,V lie on a different grid than the p points, we shall calculate the discrete derivatives ΔU/Δx(i',j') , ΔV/Δy(i',j') by mediating on the direction ortogonal to the derivative itself.
du_dx
dv_dy

Finally we can write our nu equation discretized

nu

The NS equation runs in the same staggering problem encountered in the discretization of the basin-level-equation but we should convert our nu derivative terms from the pressure grid to the velocity one





















