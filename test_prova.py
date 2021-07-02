import numpy as np
import functions as fn
from hypothesis import strategies as st
from hypothesis import given

#test functions


def custom_initialize(n):
    
    "CREATES CASUAL GRID (unitary step) WITH CUSTOM IC"
    
    if n == 0:
        
        "max velocity integral on the domain"
        
        nx = int(100* np.random.rand(1))+2
        ny = int(100* np.random.rand(1))+2    
        nz = 2
        u = np.ones(((nz,nx,ny)))
        v = np.ones(((nz,nx,ny)))
        H = np.zeros((nx+1,ny+1))
   
    elif n == 1:
        
        "max velocity divergence a point"
        np.random.seed(30)
        nx2 = int(50* np.random.rand(1))+1
        ny2 = int(50* np.random.rand(1))+1    
        nx = 2*nx2+1
        ny = 2*ny2+1
        nz = 2
        u = np.zeros(((nz,nx,ny)))
        v = np.zeros(((nz,nx,ny)))
        u[:,nx2+1,ny2]=-1
        u[:,nx2-1,ny2]=1
        v[:,nx2,ny2-1]=1
        v[:,nx2,ny2+1]=-1
        H = np.zeros((nx+1,ny+1))
    
    elif n == 2:
    
        "spike of mass in the center"
        
        np.random.seed(30)
        nx2 = int(50* np.random.rand(1))+1
        ny2 = int(50* np.random.rand(1))+1    
        nx = 2*nx2+1
        ny = 2*ny2+1
        nz = 2        
        u = np.zeros(((nz,nx,ny)))
        v = np.zeros(((nz,nx,ny)))
        H = np.zeros((nx,ny))
        H[nx2,ny2] = 2
    u[:,0,:] = 0
    v[:,0,:] = 0
    u[:,-1,:] = 0
    v[:,-1,:] = 0
    u[:,:,0] = 0
    v[:,:,0] = 0
    u[:,:,-1] = 0
    v[:,:,-1] = 0
    H[0,:] = H[1,:]
    H[-1,:]= H[-2,:]
    H[:,0] = H[:,1]
    H[:,-1] = H[:,-2]
    dt=1
    dx=1
    dy=1
    dz=1        
    return nx,ny,nz,u,v,H,dt,dx,dy,dz
        
 
def rand_initialize():

    "CREATES CASUAL GRID (unitary step) WITH CASUAL IC"
    
    np.random.seed(30)
    nx = int(100* np.random.rand(1))+2
    ny = int(100* np.random.rand(1))+2
    nz = 2
    u=np.random.rand(nz, nx, ny)
    v=np.random.rand(nz, nx, ny)
    H=2*np.random.rand(nx+1, ny+1)
    u[:,0,:] = 0
    v[:,0,:] = 0
    u[:,-1,:] = 0
    v[:,-1,:] = 0
    u[:,:,0] = 0
    v[:,:,0] = 0
    u[:,:,-1] = 0
    v[:,:,-1] = 0
    H[0,:] = H[1,:]
    H[-1,:]= H[-2,:]
    H[:,0] = H[:,1]
    H[:,-1] = H[:,-2]
    dt=1
    dx=1
    dy=1
    dz=1
    
    return nx,ny,nz,u,v,H,dt,dx,dy,dz


@given(z=st.floats(200,2000))
def test_H_time_step_mass_conservation_random_initialize(z):
    iterations = 10
    step = 0
    mH = []
    
    nx, ny, nz,u,v,H,dt,dx,dy,dz = rand_initialize()
    n = nx * ny
    sum_H0 = 0
    for i in range (1,nx):
        for j in range (1,ny):
            sum_H0 += H[i,j]
    
    
    while iterations > step :
    
        H  = fn.H_time_step(H, u, v, z, dx, dy, dt)

        sum_H1 = 0
    
        for i in range (1,nx-1):
            for j in range (1,ny-1):
                sum_H1 += H[i,j]
    

        mH.append(abs((sum_H0-sum_H1)/n))
    
        step+=1

    d_eta_max = max(mH) 
    assert 0.0000000001 > d_eta_max, "mass conservation failed"
    print("mass conservation verifed")

  

