import numpy as np
import functions as fn
from hypothesis import strategies as st
from hypothesis import given

#test functions




def custom_initialize(n):
    
    "CREATES CASUAL GRID (unitary step) WITH CUSTOM IC"
    
    if n == 0:
        
        "max wind velocity integral on the domain"
        
        nx = 10
        ny = 10
        nz = 2
        u = np.zeros(((nz,nx,ny)))
        v = np.zeros(((nz,nx,ny)))
        H = np.zeros((nx+1,ny+1))
        uw = 15*np.ones((nx,ny))
        vw = 15*np.ones((nx,ny))
        
       
   
    elif n == 1:
        
        "max wind field divergence"
        
        nx = 10
        ny = 100
        nz = 2
        u = np.zeros(((nz,nx,ny)))
        v = np.zeros(((nz,nx,ny)))
        H = np.zeros((nx+1,ny+1))
        uw = 15*np.ones((nx,ny))
        uw[5:,:]= -15
        vw = 15*np.ones((nx,ny))
        vw[:,0:]= -15
        
        
    
    elif n == 2:
    
        "max wind field rotor"
       
        nx = 10
        ny = 10 
        nz = 2
        u = np.zeros(((nz,nx,ny)))
        v = np.zeros(((nz,nx,ny)))
        H = np.zeros((nx+1,ny+1))
        uw = 15*np.ones((nx,ny))
        vw = 15*np.ones((nx,ny))
        
        uw[:,5:]= -15
        vw[5:,:]= -15
    
    elif n == 3:
        
        "CREATES CASUAL GRID (unitary step) WITH CASUAL IC"
    
        np.random.seed(30)
        nx = int(10* np.random.rand(1))+2
        ny = int(10* np.random.rand(1))+2
        nz = 2
        u=np.zeros(((nz, nx, ny)))
        v=np.zeros(((nz, nx, ny)))
        H=np.zeros(((nx+1, ny+1)))
    
        dx=100
        dy=100
        dz=100
        
        uw = 15*(np.random.rand(nx,ny))
        vw = 15*(np.random.rand(nx,ny))
    
    dt=0.9
    dx=100
    dy=100
    dz=100      
    
    return nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw

        


def test_wind_stress():
    max_dif_u = [] 
    max_dif_v = []
    for i in range(0,4):
    
        nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw = custom_initialize(i)
    
        Fx,Fy = fn.wind_stress(uw, vw)
    
        s_uw = np.sign(uw.flatten())
        s_vw = np.sign(vw.flatten())
        s_Fx = np.sign(Fx[1,:,:].flatten())
        s_Fy = np.sign(Fy[1,:,:].flatten())
    
        s_udif = abs(s_uw-s_Fx)
        s_vdif = abs(s_vw-s_Fy)
    
        max_dif_ui = max(s_udif)
        max_dif_vi = max(s_vdif)
        max_dif_u.append(max_dif_ui)
        max_dif_v.append(max_dif_vi)
    
    
    dif_u = max(max_dif_u)
    dif_v = max(max_dif_v) 

    assert dif_u == 0,""
    assert dif_v == 0,""
    
@given(nx = st.integers(2,3),ny = st.integers(2,3) )
def test_bottom_stress(nx,ny):
    u = np.random.rand(2,nx,ny)
    v = np.random.rand(2,nx,ny)
    
    Bx,By = fn.bottom_stress(u, v)
    
    s_u = np.sign(u[0,:,:].flatten())
    s_v = np.sign(v[0,:,:].flatten())
    s_Bx = np.sign(Bx[0,:,:].flatten())
    s_By = np.sign(By[0,:,:].flatten())
    
    s_udif = abs(s_u+s_Bx)
    s_vdif = abs(s_v+s_By)
    
    max_dif_u = max(s_udif)
    max_dif_v = max(s_vdif)
    

    

    assert max_dif_u == 0,""
    assert max_dif_v == 0,""


    
@given(fco = st.floats(-0.01,0.01),nu = st.floats(0,0.5),g = st.floats(9.8,9.82) )
def test_v_time_step(fco,nu,g):
    
    "verify courant condition"
    

    iterations=50
    steps=0
    dt = 0.9
    C = 0.4
    nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw = custom_initialize(3)
    z = 2*dz

    
    Fx,Fy = fn.wind_stress(uw, vw)
    
    umax = [] 
    vmax = []
    
    while iterations > steps :
        u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
        umaxi = np.max(abs(u))
        vmaxi = np.max(abs(v))
                
        umax.append(umaxi)
        vmax.append(vmaxi)
        
        steps += 1
    
    uM = np.max(umax)
    vM = np.max(vmax)
    
    uC = C*dx/dt
    vC = C*dy/dt
    
    assert uC > uM
    assert vC > vM
    


  

