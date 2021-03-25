# Closed Basin quasi-3D simulations of incompressible unsteady flow in a closed basin with idrostatic approx

PHISICAL FORMULATION

In order to simulate the behavior of a fluid in a close basin we should start from the navier stokes equations:

# ∂u∂t+u∂u∂x+v∂u∂y=−1ρ∂p∂x+ν(∂2u∂x2+∂2u∂y2)+Fx
# ∂v∂t+u∂v∂x+v∂v∂y=−1ρ∂p∂y+ν(∂2v∂x2+∂2v∂y2)+Fy
# ∂w∂t+u∂v∂x+v∂v∂y=−1ρ∂p∂y+ν(∂2v∂x2+∂2v∂y2)+Fz

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

As boundary condition for a closed basin seems reasonable a no-slip condition and a null pressure gradient
u(boundary)=0
v(boundary)=0
∂p∂xi(boundary)=0


the system is closed and ,given an initial condition, we are ready to numerically solve it.

NUMERICAL FORMULATION

DEFINITION OF THE GRID

The grids used are the easiest possible: uniform 3-D  solidal with our x,y,z axes.
The only complication is the using of two grids, one for the velocity and one staggered grid for the pressure, with each cell of the second one centered in the intesections of the cells of the first one.

For general purpose and definition of our algorithms, a Staggered CFD Grid shall be used. It ensures the presence of very real non-zero pressure gradient across the nodes in any condition, even in the case of a checkered grid. The staggering also ensures realistic behaviors of the descretized momentum equations for spatially oscillating pressures. Also, the direction of the velocity vectors are exact.

DISCRETIZATION OF THE EQUATIONS

  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "$$\\frac{\\partial u}{\\partial t}+u\\frac{\\partial u}{\\partial x}+v\\frac{\\partial u}{\\partial y}=-\\frac{1}{\\rho}\\frac{\\partial p}{\\partial x}+\\nu\\left(\\frac{\\partial^2 u}{\\partial x^2}+\\frac{\\partial^2 u}{\\partial y^2}\\right)+F$$\n",
    "\n",
    "$$\\frac{\\partial v}{\\partial t}+u\\frac{\\partial v}{\\partial x}+v\\frac{\\partial v}{\\partial y}=-\\frac{1}{\\rho}\\frac{\\partial p}{\\partial y}+\\nu\\left(\\frac{\\partial^2 v}{\\partial x^2}+\\frac{\\partial^2 v}{\\partial y^2}\\right)$$\n",
    "\n",
    "$$\\frac{\\partial^2 p}{\\partial x^2}+\\frac{\\partial^2 p}{\\partial y^2}=-\\rho\\left(\\frac{\\partial u}{\\partial x}\\frac{\\partial u}{\\partial x}+2\\frac{\\partial u}{\\partial y}\\frac{\\partial v}{\\partial x}+\\frac{\\partial v}{\\partial y}\\frac{\\partial v}{\\partial y}\\right)\n",
    "$$"
   ]
  },














