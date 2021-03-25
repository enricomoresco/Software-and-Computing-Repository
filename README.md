# Closed Basin quasi-3D simulations of incompressible unsteady flow in a closed basin with idrostatic approx

PHISICAL FORMULATION

In order to simulate the behavior of a fluid in a close basin we should start from the navier stokes equations:

∂u∂t+u∂u∂x+v∂u∂y=−1ρ∂p∂x+ν(∂2u∂x2+∂2u∂y2)+Fx
∂v∂t+u∂v∂x+v∂v∂y=−1ρ∂p∂y+ν(∂2v∂x2+∂2v∂y2)+Fy
∂w∂t+u∂v∂x+v∂v∂y=−1ρ∂p∂y+ν(∂2v∂x2+∂2v∂y2)+Fz

Fx, Fy, Fz in general are the resultant of any external force acting on our system:
in our case we have the gravity force acting along the z-axis 
Fz=-ρg

and the coriolis inertial force acting along x and y axes as result of earth rotation (the centrifugal force is neglectable)

Fx=fv
Fy=-fu

f is the coriolis parameter, is about 10E-4 rad/s mutiplied by the sine of the latitude.
For a wide shallow basin we can simplify the third equation with the hidrostatic approximation:
p=-ρgh∂n∂x

Considerin a uniform water density we can write the spatial derivatives of p, for example ∂p∂x as
∂p∂x=-ρg∂n∂x

with n the elevation of the surface of the basin

Then we consider incompressibility

∂u∂x+∂v∂y+∂w∂z=0

Integrating on all the water column we obtain
-(∂U∂x+∂V∂y)=∂n∂t

with U,V integrals of the velocity over all the water column 

The system of differential equations now looks like this

∂u∂t+u∂u∂x+v∂u∂y=−1ρ∂p∂x+ν(∂2u∂x2+∂2u∂y2)+fv
∂v∂t+u∂v∂x+v∂v∂y=−1ρ∂p∂y+ν(∂2v∂x2+∂2v∂y2)-fu
∂p∂xi=-ρg∂n∂xi
-(∂U∂x+∂V∂y)=∂n∂t

the system is closed and we are ready to numerically solve it.

NUMERICAL FORMULATION















