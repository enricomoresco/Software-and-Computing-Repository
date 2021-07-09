import numpy
import copy as cp

#calculate bottom stress

def bottom_stress(u, v):
    """this function calculates the bottom stress from the velocity field using
    a default bottom rugosity constant.
    
    input: u(zonal current velocity), v(meridional current velocity)
    output: Bx(zonal bottom stress), By(meridional bottom stress)"""
    
    nx = len(u[0,:,0])
    ny = len(u[0,0,:])
    nz = 2
    Bx = numpy.zeros(((nz,nx,ny)))
    By = numpy.zeros(((nz,nx,ny)))
    k = 0.01
    Bx[0,:,:]= -k*u[0,:,:]*numpy.sqrt((u[0,:,:]**2)+(v[0,:,:]**2))
    By[0,:,:]= -k*v[0,:,:]*numpy.sqrt((u[0,:,:]**2)+(v[0,:,:]**2))
    return Bx, By

#calculate wind stress

def wind_stress(uw, vw):
    """this function calculates the wind stress from the velocity field using
    a default air-water density-rate constant.
    
    input: u(zonal wind velocity), v(meridional wind velocity)
    output: Bx(zonal wind stress), By(meridional wind stress)"""
    
    nx = len(uw[:,0])
    ny = len(uw[0,:])
    nz = 2    
    Fx = numpy.zeros(((nz,nx,ny)))
    Fy = numpy.zeros(((nz,nx,ny)))
    k = 0.001
    Fx[1,:,:]= k*uw[:,:]*numpy.sqrt((uw[:,:]**2)+(vw[:,:]**2))
    Fy[1,:,:]= k*vw[:,:]*numpy.sqrt((uw[:,:]**2)+(vw[:,:]**2))
    return Fx, Fy


#x-derivative of 2-variables function

def Dexb(f, dx):
    nx = len(f[:,0])
    ny = len(f[0,:])
    
    f_1= numpy.zeros((nx,ny))
    f_1[:-1,:]=(f[1:,:]-f[:-1,:])/dx
    return f_1

#y-derivative of 2-variables function

def Deyb(f, dy):
    nx = len(f[:,0])
    ny = len(f[0,:])
    f_1 =  numpy.zeros((nx,ny))
    f_1[:,:-1] = (f[:,1:]-f[:,:-1])/dy
    return f_1

#x-derivative of 3-variables function

def Dex(f, dx):
    nx = len(f[0,:,0])
    ny = len(f[0,0,:])
    nz = 2    
    f_1 = numpy.zeros(((nz,nx,ny)))
    f_1[:,:-1,:] = (f[:,1:,:]-f[:,:-1,:])/dx
    return f_1

#y-derivative of 3-variables function

def Dey(f, dy):
    nx = len(f[0,:,0])
    ny = len(f[0,0,:])
    nz = 2  
    f_1 = numpy.zeros(((nz,nx,ny)))
    f_1[:,:,:-1] = (f[:,:,1:]-f[:,:,:-1])/dy
    return f_1

#second x-derivative of 3-variables function
def Dex2(f, dx):
    nx = len(f[0,:,0])
    ny = len(f[0,0,:])
    nz = 2  
    f_2 = numpy.zeros(((nz,nx,ny)))
    f_2[:,1:-1,:] = (f[:,2:nx,:]+f[:,0:nx-2,:]-(2*f[:,1:-1,:]))/dx**2
    return f_2

#second y-derivative of 3-variables function

def Dey2(f, dy):
    nx = len(f[0,:,0])
    ny = len(f[0,0,:])
    nz = 2  
    f_2 = numpy.zeros(((nz,nx,ny)))
    f_2[:,:,1:-1] = (f[:,:,2:ny]+f[:,:,0:ny-2]-(2*f[:,:,1:-1]))/dy**2
    return f_2

#second x-derivative of 2-variables function
def Dex2b(f,dx):
    nx = len(f[:,0])
    ny = len(f[0,:])
    f_2 = numpy.zeros((nx,ny))
    f_2[1:-1,:] = (f[2:,:]+f[0:-2,:]-(2*f[1:-1,:]))/dx**2
    return f_2

#second y-derivative of 2-variables function

def Dey2b(f, nx, ny, dy):
    nx = len(f[:,0])
    ny = len(f[0,:])
    f_2 = numpy.zeros((nx,ny))
    f_2[:,1:-1] = (f[:,2:]+f[:,0:-2]-(2*f[:,1:-1]))/dy**2
    return f_2

#Advection terms

"""the following functions are the four advection terms needed for the
implementation of the NS-equations,
input : 
        u(zonal current velocity), v(meridional current velocity ),
        dx and dy (grid steps)

output : udxu,udxv,vdyu,vdyv (the four advection terms)
"""

def udexu(u,dx):
    un = cp.deepcopy(u)  
    Dexun = Dex(un,dx)
    nx = len(u[0,:,0])
    ny = len(u[0,0,:])
    nz = 2 
    udxu = numpy.zeros(((nz,nx,ny)))
    udxu[:,1:-1,1:-1] = un[:,1:-1,1:-1]*(Dexun[:,:-2,1:-1]+Dexun[:,1:-1,1:-1])/2
    
    return udxu


def udexv(u,v,dx):
    un = cp.deepcopy(u)
    vn = cp.deepcopy(v)
    Dexvn = Dex(vn,dx)
    nx = len(u[0,:,0])
    ny = len(u[0,0,:])
    nz = 2 
    udxv = numpy.zeros(((nz,nx,ny)))
    udxv[:,1:-1,1:-1] = un[:,1:-1,1:-1]*(Dexvn[:,:-2,1:-1]+Dexvn[:,1:-1,1:-1])/2
    return udxv

def vdeyu(u,v,dy):
    un = cp.deepcopy(u)
    vn = cp.deepcopy(v)
    Deyun = Dey(un,dy)  
    nx = len(u[0,:,0])
    ny = len(u[0,0,:])
    nz = 2 
    vdyu = numpy.zeros(((nz,nx,ny)))
    vdyu[:,1:-1,1:-1] = vn[:,1:-1,1:-1]*(Deyun[:,1:-1,:-2]+Deyun[:,1:-1,1:-1])/2
    return vdyu

def vdeyv(v,dy):
    nx = len(v[0,:,0])
    ny = len(v[0,0,:])
    nz = 2 
    vdyv = numpy.zeros(((nz,nx,ny)))
    vn = cp.deepcopy(v)
    Deyvn = Dey(vn,dy)
    vdyv[:,1:-1,1:-1] = vn[:,1:-1,1:-1]*(Deyvn[:,1:-1,:-2]+Deyvn[:,1:-1,1:-1])/2
    return vdyv
   
#Define how the elevation evolves through time 



def H_time_step(H,u,v,z,dx,dy,dt):
    
    """this function implements the evolution of eta in the time step
    from the vertical integral velocities U and V , using incompressibility and continuity equation
    input: H(elevation of the surface), u and v (current horizontal velocities), dx and dy (grid steps), dt(t step)
    output:H(elevation of the surface)
    
    for futher deepening is provided "Phisical_and_Numerical_formulation"
    """
    nx = len(u[0,:,0])
    ny = len(u[0,0,:])
    nz = 2 
    Hn = cp.deepcopy(H)
    U= numpy.zeros((nx+1,ny+1))
    V= numpy.zeros((nx+1,ny+1))
    U[1:,1:]=(sum(u[:,:,:]))*(z+Hn[:-1,:-1])/nz
    V[1:,1:]=(sum(v[:,:,:]))*(z+Hn[:-1,:-1])/nz
   
    DexbU = Dexb(U,dx)
    DeybV = Deyb(V,dy)
    H[1:-1,1:-1]=Hn[1:-1,1:-1]-dt*((DexbU[1:-1,1:-1]+DexbU[1:-1,2:])/2+(DeybV[1:-1,1:-1]+DeybV[2:,1:-1])/2)
    #BC gradiente di pressione nullo al bordo lungo la perpendicolare
    H[:,0] = H[:,1]
    H[:,ny]=H[:,ny-1]
    H[0,:]  = H[1,:]
    H[nx,:] = H[nx-1,:]

    return H

#Define how the velocity evolves through time 

def vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu):
    
    """this function implements the evolution of all the variables (u,v and H) 
    in the time step:
    for H simply uses "H_time_step" as implemented before
    for u, v uses a simplification of NS equation
    
    input: u and v (horizontal current velocities), H(elevation of the surface),
    Fx and Fy (horizontal wind stress), dx and dy (grid steps), dt(t step), 
    g(gravity acceleration), fco(coriolis factor), nu(eddy viscosity)
    
    output: u and v (horizontal current velocities), H(elevation of the surface)
    udiff,vdiff,Hdiff (elements for the implementation of stationary solution research)
    
    for futher deepening is provided "Phisical_and_Numerical_formulation"
    """
    nx = len(u[0,:,0])
    ny = len(u[0,0,:])
    nz = 2 
    Hn = H.copy()
    H = H_time_step(H,u,v,z,dx,dy,dt)
    
    Bx,By = bottom_stress(u, v)
    
    cox = numpy.zeros(((nz,nx,ny)))
    coy = numpy.zeros(((nz,nx,ny)))
    dexP = numpy.zeros((nx,ny))
    deyP = numpy.zeros((nx,ny))

    disu = numpy.zeros(((nz,nx,ny)))
    disv = numpy.zeros(((nz,nx,ny)))
    Dez2un = numpy.zeros(((nz,nx,ny)))
    Dez2vn = numpy.zeros(((nz,nx,ny)))
    
    un = u.copy()
    vn = v.copy()

    Dez2un[0,:,:]=-(un[0,:,:]-un[1,:,:])/(dz**2)
    Dez2un[1,:,:]=-Dez2un[0,:,:]
    Dez2vn[0,:,:]=-(vn[0,:,:]-vn[1,:,:])/(dz**2)
    Dez2vn[1,:,:]=-Dez2vn[0,:,:]
     
    
    cox[:,:,:] = fco*vn[:,:,:]
    coy[:,:,:] = -fco*un[:,:,:]
    udxu = udexu(u, dx)
    udxv = udexv(u,v, dx)
    vdyu = vdeyu(u,v, dy)
    vdyv = vdeyv(v, dy)
    dexP[:,:] = g/2 * (Dexb(H,dx)[:-1,:-1]+Dexb(H,dx)[:-1,1:])
    deyP[:,:] = g/2 * (Deyb(H,dy)[:-1,:-1]+Deyb(H,dy)[1:,:-1])
    disuh = nu * (Dex2(un,dx) + Dey2(un,dy))
    disvh = nu * (Dex2(vn,dx) + Dey2(vn,dy))
    disu[:,:,:] = disuh[:,:,:] + Dez2un[:,:,:]
    disv[:,:,:] = disvh[:,:,:] + Dez2vn[:,:,:]
    
    u[:,1:-1,1:-1] = (un[:,1:-1,1:-1] - dexP[1:-1,1:-1]-udxu[:,1:-1,1:-1]-vdyu[:,1:-1,1:-1]+disu[:,1:-1,1:-1]+cox[:,1:-1,1:-1]+Fx[:,1:-1,1:-1]+Bx[:,1:-1,1:-1])*dt
    v[:,1:-1,1:-1] = (vn[:,1:-1,1:-1] - deyP[1:-1,1:-1]-udxv[:,1:-1,1:-1]-vdyv[:,1:-1,1:-1]+disv[:,1:-1,1:-1]+coy[:,1:-1,1:-1]+Fy[:,1:-1,1:-1]+By[:,1:-1,1:-1])*dt

    du4 = (u-un)**4
    dv4 = (v-vn)**4
    dH2 = (H-Hn)**2
   
    u4 = u**4
    v4 = v**4
    H2 = H**2
    g2 = g**2

    udiff = numpy.sum(du4)/(numpy.sum(u4)+numpy.sum(v4)+g2*numpy.sum(H2))
    vdiff = numpy.sum(dv4)/(numpy.sum(u4)+numpy.sum(v4)+g2*numpy.sum(H2))
    Hdiff = numpy.sum(dH2)/(numpy.sum(H2)+numpy.sum(u4)/g2+numpy.sum(v4)/100)
    
    return u,v,H,udiff,vdiff,Hdiff

