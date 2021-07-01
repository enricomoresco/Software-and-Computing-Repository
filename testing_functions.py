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
        
        "max velocity divergence in the domain"
        
        nx2 = int(50* np.random.rand(1))+1
        ny2 = int(50* np.random.rand(1))+1    
        nx = 2*nx2+1
        ny = 2*nx2+1
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
        
        nx2 = int(50* np.random.rand(1))+1
        ny2 = int(50* np.random.rand(1))+1    
        nx = 2*nx2+1
        ny = 2*nx2+1
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
    
    while iterations > step :
    
        nx, ny, nz,u,v,H,dt,dx,dy,dz = rand_initialize()
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

@given(z=st.floats(200,2000))    
def test_H_timestep_mass_conservation_custom_initialize(z):

    mH = []
    
    for i in range (0,2):
    
        nx, ny, nz,u,v,H,dt,dx,dy,dz = custom_initialize(i)
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
    
       

    d_eta_max = max(mH) 
    assert 0.0000000001 > d_eta_max, "mass conservation failed"
    print("mass conservation verifed")

#to verify the momentum conservation in the advection parallel terms:


def test_udxu_x_momentum_conservation_random_initialize(iterations=10,step=0):
    
    mx = []    
    
    while iterations > step :
        


        nx, ny, nz,u,v,H,dt,dx,dy,dz = rand_initialize()
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
    


def test_vdyvy_momentum_conservation_random_initialize(iterations=10,step=0):
    
    my = []
    
    
    while iterations > step :
        


        nx, ny, nz,u,v,H,dt,dx,dy,dz = rand_initialize()
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
    
    

def test_vdyv_y_momentum_conservation_custom():
    
    my = []
    
    
    for i in range (0,2):
        
        nx, ny, nz,u,v,H,dt,dx,dy,dz = custom_initialize(i)
        n = nx * ny
        vdyv = fn.vdeyv(v, dy)
        sum_adv_v=0
        
        for i in range (1,nx-1):
            for j in range (1,ny-1):
                for k in range (0,2):
                    
                    sum_adv_v+=vdyv[k,i,j]

        my.append(sum_adv_v/n)
        
    
    d_v_max = max(my)
    assert 0.0000000001 > d_v_max, "y-momentum conservation failed"
    
    print("y-momentum conservation verifed")



@given(z=st.floats(200,2000),fco=st.floats(-0.0001,0.0001),nu =st.floats(0.0001,0.01),g=st.floats(9.80,9.82)) 
def test_vel_time_step_u_BC_evolution_random(z,fco,nu,g):
    iterations=10
    step=0
    while iterations > step :
    
        nx, ny, nz,u,v,H,dt,dx,dy,dz = rand_initialize()
    
        Fx = np.zeros(((nz,nx,ny)))
        Fy = np.zeros(((nz,nx,ny)))
        
        u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
       
        u_abs_x_0 = abs(u[:,0,:])
        u_abs_x_n = abs(u[:,-1,:])
        u_abs_y_0 = abs(u[:,:,0])
        u_abs_y_n = abs(u[:,:,-1])
        
        uabs_x_0 = u_abs_x_0.flatten()
        uabs_x_n = u_abs_x_n.flatten()
        uabs_y_0 = u_abs_y_0.flatten()
        uabs_y_n = u_abs_y_n.flatten()
               
        umax_x_0 = max(uabs_x_0)
        umax_x_n = max(uabs_x_n)
        umax_y_0 = max(uabs_y_0)
        umax_y_n = max(uabs_y_n)
        
        assert umax_x_0 == 0, "u BC evolution failed"
        assert umax_x_n == 0, "u BC evolution failed"
        assert umax_y_0 == 0, "u BC evolution failed"
        assert umax_y_n == 0, "u BC evolution failed"
        
        step += 1
    
    print("u BC evolution verified")
        


@given(z=st.floats(200,2000),fco=st.floats(-0.0001,0.0001),nu =st.floats(0.0001,0.01),g=st.floats(9.80,9.82)) 
def test_vel_time_step_u_BC_evolution_custom_initialize(z,fco,nu,g):
    for i in range (0,2):
    
        nx, ny, nz,u,v,H,dt,dx,dy,dz = custom_initialize(i)
    
        Fx = np.zeros(((nz,nx,ny)))
        Fy = np.zeros(((nz,nx,ny)))
        
        u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
       
        u_abs_x_0 = abs(u[:,0,:])
        u_abs_x_n = abs(u[:,-1,:])
        u_abs_y_0 = abs(u[:,:,0])
        u_abs_y_n = abs(u[:,:,-1])
        
        uabs_x_0 = u_abs_x_0.flatten()
        uabs_x_n = u_abs_x_n.flatten()
        uabs_y_0 = u_abs_y_0.flatten()
        uabs_y_n = u_abs_y_n.flatten()
               
        umax_x_0 = max(uabs_x_0)
        umax_x_n = max(uabs_x_n)
        umax_y_0 = max(uabs_y_0)
        umax_y_n = max(uabs_y_n)
        
        assert umax_x_0 == 0, "u BC evolution failed"
        assert umax_x_n == 0, "u BC evolution failed"
        assert umax_y_0 == 0, "u BC evolution failed"
        assert umax_y_n == 0, "u BC evolution failed"
        
        print("u BC evolution verified")
    


@given(z=st.floats(200,2000),fco=st.floats(-0.0001,0.0001),nu =st.floats(0.0001,0.01),g=st.floats(9.80,9.82)) 
def test_vel_time_step_v_BC_evolution_random_initialize(z,fco,nu,g):
    iterations=10
    step=0
    while iterations > step :
    
        nx, ny, nz,u,v,H,dt,dx,dy,dz = rand_initialize()
    
        Fx = np.zeros(((nz,nx,ny)))
        Fy = np.zeros(((nz,nx,ny)))
        
        u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
       
        v_abs_x_0 = abs(v[:,0,:])
        v_abs_x_n = abs(v[:,-1,:])
        v_abs_y_0 = abs(v[:,:,0])
        v_abs_y_n = abs(v[:,:,-1])
        
        vabs_x_0 = v_abs_x_0.flatten()
        vabs_x_n = v_abs_x_n.flatten()
        vabs_y_0 = v_abs_y_0.flatten()
        vabs_y_n = v_abs_y_n.flatten()
        
        
        vmax_x_0 = max(vabs_x_0)
        vmax_x_n = max(vabs_x_n)
        vmax_y_0 = max(vabs_y_0)
        vmax_y_n = max(vabs_y_n)
        
        assert vmax_x_0 == 0, "v BC evolution failed"
        assert vmax_x_n == 0, "v BC evolution failed"
        assert vmax_y_0 == 0, "v BC evolution failed"
        assert vmax_y_n == 0, "v BC evolution failed"

        step += 1
        
@given(z=st.floats(200,2000),fco=st.floats(-0.0001,0.0001),nu =st.floats(0.0001,0.01),g=st.floats(9.80,9.82)) 
def test_vel_time_step_v_BC_evolution_custom_initialize(z,fco,nu,g):
    for i in range (0,2):
    
        nx, ny, nz,u,v,H,dt,dx,dy,dz = custom_initialize(i)
    
        Fx = np.zeros(((nz,nx,ny)))
        Fy = np.zeros(((nz,nx,ny)))
        
        u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt,g,fco,nu)
       
        v_abs_x_0 = abs(v[:,0,:])
        v_abs_x_n = abs(v[:,-1,:])
        v_abs_y_0 = abs(v[:,:,0])
        v_abs_y_n = abs(v[:,:,-1])
        
        vabs_x_0 = v_abs_x_0.flatten()
        vabs_x_n = v_abs_x_n.flatten()
        vabs_y_0 = v_abs_y_0.flatten()
        vabs_y_n = v_abs_y_n.flatten()
        
        
        vmax_x_0 = max(vabs_x_0)
        vmax_x_n = max(vabs_x_n)
        vmax_y_0 = max(vabs_y_0)
        vmax_y_n = max(vabs_y_n)
        
        assert vmax_x_0 == 0, "v BC evolution failed"
        assert vmax_x_n == 0, "v BC evolution failed"
        assert vmax_y_0 == 0, "v BC evolution failed"
        assert vmax_y_n == 0, "v BC evolution failed"

        
    
    print("v BC evolution verified")
    
    print("v BC evolution verified")
    
@given(z=st.floats(200,2000),fco=st.floats(-0.0001,0.0001),nu =st.floats(0.0001,0.01),g=st.floats(9.80,9.82)) 
def test_vel_time_step_H_BC_evolution_random_initialize(z,fco,nu,g):
    iterations=10
    step=0
    while iterations > step :
    
        nx, ny, nz,u,v,H,dt,dx,dy,dz = rand_initialize()
    
        Fx = np.zeros(((nz,nx,ny)))
        Fy = np.zeros(((nz,nx,ny)))
        
        u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz,dt,g,fco,nu)
       
        Hdiff_x_0 = abs(H[0,1:-1] - H[1,1:-1])
        Hdiff_x_n = abs(H[-1,1:-1]- H[-2,1:-1])
        Hdiff_y_0 = abs(H[1:-1,0] - H[1:-1,1])
        Hdiff_y_n = abs(H[1:-1,-1] - H[1:-1,-2])
        
        
        Hdmax_x_0 = max(Hdiff_x_0)
        Hdmax_x_n = max(Hdiff_x_n)
        Hdmax_y_0 = max(Hdiff_y_0)
        Hdmax_y_n = max(Hdiff_y_n)
        
        assert Hdmax_x_0 == 0, "H BC evolution failed"
        assert Hdmax_x_n == 0, "H BC evolution failed"
        assert Hdmax_y_0 == 0, "H BC evolution failed"
        assert Hdmax_y_n == 0, "H BC evolution failed"

        step += 1
    
    print("H BC evolution verified")
    
@given(z=st.floats(200,2000),fco=st.floats(-0.0001,0.0001),nu =st.floats(0.0001,0.01),g=st.floats(9.80,9.82)) 
def test_vel_time_step_H_BC_evolution_custom_initialize(z,fco,nu,g):
    for i in range (0,2):
    
    
        nx, ny, nz,u,v,H,dt,dx,dy,dz = custom_initialize(i)
    
        Fx = np.zeros(((nz,nx,ny)))
        Fy = np.zeros(((nz,nx,ny)))
        
        u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz,dt,g,fco,nu)
       
        Hdiff_x_0 = abs(H[0,1:-1] - H[1,1:-1])
        Hdiff_x_n = abs(H[-1,1:-1]- H[-2,1:-1])
        Hdiff_y_0 = abs(H[1:-1,0] - H[1:-1,1])
        Hdiff_y_n = abs(H[1:-1,-1] - H[1:-1,-2])
        
        
        Hdmax_x_0 = max(Hdiff_x_0)
        Hdmax_x_n = max(Hdiff_x_n)
        Hdmax_y_0 = max(Hdiff_y_0)
        Hdmax_y_n = max(Hdiff_y_n)
        
        assert Hdmax_x_0 == 0, "H BC evolution failed"
        assert Hdmax_x_n == 0, "H BC evolution failed"
        assert Hdmax_y_0 == 0, "H BC evolution failed"
        assert Hdmax_y_n == 0, "H BC evolution failed"

        
    
    print("H BC evolution verified")
    



#to verify the momentum conservation in the advection parallel terms:


def test_udxu_x_momentum_conservation_custom_initialize():
    
    mx = []    
    
    for i in range (0,2):
        


        nx, ny, nz,u,v,H,dt,dx,dy,dz = custom_initialize(i)
        n = nx * ny

        udxu = fn.udexu(u, dx)
        sum_adv_u=0
        
        for i in range (1,nx-1):
            for j in range (1,ny-1):
                for k in range (0,2):
                    
                    sum_adv_u+=udxu[k,i,j]



        mx.append(sum_adv_u/n)


        
    
    
    d_u_max = max(mx)


        
    assert 0.0000000001 > d_u_max, "x-momentum conservation failed"

    
    print("x-momentum conservation verifed")
    





    

    

    
                
        


