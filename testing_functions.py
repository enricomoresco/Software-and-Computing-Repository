import numpy as np
import functions as fn
from hypothesis import strategies as st
from hypothesis import given
from hypothesis import settings


def custom_initialize(n):
    
    "CREATES CASUAL GRID (unitary step) WITH CUSTOM IC"
    
    if n == 0:
        
        "max wind velocity integral on the domain"
        
        nx = 10
        ny = 10
        nz = 2

        uw = 15*np.ones((nx,ny))
        vw = 15*np.ones((nx,ny))
        
       
   
    elif n == 1:
        
        "max wind field divergence"
        
        nx = 10
        ny = 10
        nz = 2

        uw = 15*np.ones((nx,ny))
        uw[5:,:]= -15
        vw = 15*np.ones((nx,ny))
        vw[:,0:]= -15
        
        
    
    elif n == 2:
    
        "max wind field rotor"
       
        nx = 10
        ny = 10 
        nz = 2

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

    
        dx=100
        dy=100
        dz=100
        
        uw = 15*(np.random.rand(nx,ny))
        vw = 15*(np.random.rand(nx,ny))
   
    u=np.zeros(((nz, nx, ny)))
    v=np.zeros(((nz, nx, ny)))
    H=np.zeros(((nx+1, ny+1)))
    dt=0.9
    dx=100
    dy=100
    dz=100      
    
    return nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw





#test functions




  
def test_H_time_step():
    "this function tests the mass conservation during the iterations"
    iterations = 50
    step = 0
    mH = []
    for n in range (3):
        nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw = custom_initialize(n)
    
        while iterations > step :
    
        
            z = 2* dz
            n = nx * ny
            sum_H0 = 0

            for i in range (1,nx):
                for j in range (1,ny):
                    sum_H0 += H[i,j]

            H = fn.H_time_step(H,u,v,z,dx,dy,dt) 
            
            sum_H1 = 0
    
            for i in range (1,nx):
                for j in range (1,ny):
                    sum_H1 += H[i,j]
    

            mH.append(abs((sum_H0-sum_H1)/n))
    
            step+=1

    d_eta_max = max(mH) 
    assert 0.0000000001 > d_eta_max, "mass conservation failed"
    print("mass conservation verifed")



@settings(max_examples=50,deadline =500)
@given(fco = st.floats(-0.01,0.01),nu = st.floats(0,0.5),g = st.floats(9.8,9.82))
def test_udexu(fco,nu,g):
    "this funciton tests the conservation of the x-momentum trough the advecton effect"
    
    iterations = 50
    step = 0
    
    
    mx = []    
    for n in range (3):
        nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw = custom_initialize(n)
        Fx,Fy = fn.wind_stress(uw, vw)
        z = 2*dz
        u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
    
        while iterations > step :
        
            u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
        
            n = nx * ny

            udxu = fn.udexu(u, dx)
            sum_adv_u=0
        
            for i in range (1,nx-1):
                for j in range (1,ny-1):
                    for k in range (0,2):
                    
                        sum_adv_u+=udxu[k,i,j]



            mx.append(sum_adv_u/n)


            step+=1
    
    
    d_u_max = max(mx)


        
    assert 0.0000000001 > d_u_max, "x-momentum conservation failed"

    
    print("x-momentum conservation verifed")
    


@settings(max_examples=50,deadline =500) 
@given(fco = st.floats(-0.01,0.01),nu = st.floats(0,0.5),g = st.floats(9.8,9.82))
def test_vdeyv(fco,nu,g):
    "this funciton tests the conservation of the y-momentum trough the advecton effect"
    
    iterations = 50
    step = 0
    
    my = []
    
    for n in range (3):
        nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw = custom_initialize(n)
        Fx,Fy = fn.wind_stress(uw, vw)
        z = 2*dz
        u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
    
    
        while iterations > step :
            
            u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
            
            n = nx * ny

            vdyv = fn.vdeyv(v, dy)
        
            sum_adv_v=0
        
            for i in range (1,nx-1):
                for j in range (1,ny-1):
                    for k in range (0,2):
                    
                        sum_adv_v+=vdyv[k,i,j]

            my.append(sum_adv_v/n)
            step+=1
    
    d_v_max = max(my)
    assert 0.0000000001 > d_v_max, "y-momentum conservation failed"
    
    print("y-momentum conservation verifed")
    
    


"test of the evolution of the boundary conditions for the function vel_time_step"

@settings(max_examples=50,deadline =500)
@given(fco=st.floats(-0.0001,0.0001),nu =st.floats(0.0001,0.01),g=st.floats(9.80,9.82)) 
def test_vel_time_step_BC(fco,nu,g):
    """this fuction test the persistence of Boundary Conditions physically acceptable
    (no mass or energy trasmission trough the edge)"""
    iterations=50
    step=0
    
    for n in range (3):
        nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw = custom_initialize(n)
        Fx,Fy = fn.wind_stress(uw, vw)
        z = 2*dz
        while iterations > step :

        
            u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
       

            umax_x_0 = np.max(abs(u[:,0,:]))
            umax_x_n = np.max(abs(u[:,-1,:]))
            umax_y_0 = np.max(abs(u[:,:,0]))
            umax_y_n = np.max(abs(u[:,:,-1]))
            vmax_x_0 = np.max(abs(v[:,0,:]))
            vmax_x_n = np.max(abs(v[:,-1,:]))
            vmax_y_0 = np.max(abs(v[:,:,0]))
            vmax_y_n = np.max(abs(v[:,:,-1]))
        
            assert umax_x_0 == 0, "u BC evolution failed"
            assert umax_x_n == 0, "u BC evolution failed"
            assert umax_y_0 == 0, "u BC evolution failed"
            assert umax_y_n == 0, "u BC evolution failed"
            assert vmax_x_0 == 0, "u BC evolution failed"
            assert vmax_x_n == 0, "u BC evolution failed"
            assert vmax_y_0 == 0, "u BC evolution failed"
            assert vmax_y_n == 0, "u BC evolution failed"
            step += 1
    
    print("BC evolution verified")
        

def test_wind_stress():
    
    """"test wind stress parallel direction respect to the wind velocity field
    (wind positive work)"""
     
    max_dif_u = [] 
    max_dif_v = []
    
    for n in range (3):
        nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw = custom_initialize(n)
    
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

@settings(max_examples=50,deadline =500)
@given(fco=st.floats(-0.0001,0.0001),nu =st.floats(0.0001,0.01),g=st.floats(9.80,9.82)) 
def test_bottom_stress(fco,nu,g):
    
    """test bottom stress opposite direction respect to the current velocity field
    (dissipation)"""
    
    
    iterations=50
    step=0
    max_dif_u = [] 
    max_dif_v = []
    
    for n in range (3):
        nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw = custom_initialize(n)
        Fx,Fy = fn.wind_stress(uw, vw)
        z = 2*dz

        while iterations > step :
    
            nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw = custom_initialize(n)
    
            Bx,By = fn.bottom_stress(u, v)
    
            s_u = np.sign(u[0,:,:].flatten())
            s_v = np.sign(v[0,:,:].flatten())
            s_Bx = np.sign(Bx[0,:,:].flatten())
            s_By = np.sign(By[0,:,:].flatten())
    
            s_udif = abs(s_u+s_Bx)
            s_vdif = abs(s_v+s_By)
    
            max_dif_ui = max(s_udif)
            max_dif_vi = max(s_vdif)
        
            max_dif_u.append(max_dif_ui)
            max_dif_v.append(max_dif_vi)
        
            u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
    
            step += 1
    
    dif_u = max(max_dif_u)
    dif_v = max(max_dif_v) 

    assert dif_u == 0,""
    assert dif_v == 0,""
        
        
        
@settings(max_examples=50,deadline =500)
@given(fco = st.floats(-0.01,0.01),nu = st.floats(0,0.5),g = st.floats(9.8,9.82))
def test_v_time_step_courant(fco,nu,g):
    
    """verify courant stability condition:
        time step * courant constant >  grid step / velocity"""
    

    iterations=50
    steps=0
    dt = 0.9
    C = 0.4
    umax = [] 
    vmax = []
    
    for n in range (3):
        nx,ny,nz,u,v,H,dt,dx,dy,dz,uw,vw = custom_initialize(n)
        z = 2*dz

    
        Fx,Fy = fn.wind_stress(uw, vw)
    

    
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







    

    

    
                
        


