import numpy as np
import functions as fn


#test functions

def custom_inizialize(n):
    if n == 0:
        #max velocity module integral on the domain
        nx = int(100* np.random.rand(1))+2
        ny = int(100* np.random.rand(1))+2    
        nz = 2
        u = np.ones(((nz,nx,ny)))
        v = np.ones(((nz,nx,ny)))
        H = np.zeroes((nx+1,ny+1))
    elif n == 1:
        #max velocity divergence in the domain
        nx2 = int(50* np.random.rand(1))+1
        ny2 = int(50* np.random.rand(1))+1    
        nx = 2*nx2+1
        ny = 2*nx2+1
        nz = 2
        u = np.zeroes(((nz,nx,ny)))
        v = np.zeroes(((nz,nx,ny)))
        u[:,nx2+1,ny2]=-1
        u[:,nx2-1,ny2]=1
        v[:,nx2,ny2-1]=1
        v[:,nx2,ny2+1]=-1
        H = np.zeroes((nx+1,ny+1))
    elif n == 2:
        #spike of mass in the center
        nx2 = int(50* np.random.rand(1))+1
        ny2 = int(50* np.random.rand(1))+1    
        nx = 2*nx2+1
        ny = 2*nx2+1
        nz = 2        
        u = np.zeroes(((nz,nx,ny)))
        v = np.zeroes(((nz,nx,ny)))
        H = np.zeroes((nx,ny))
        H[nx2,ny2] = 2
        
 

        
        
    

#CREATES CASUAL GRID (unitary step)vWITH CASUAL IC
def rand_initialize():
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
   

def mass_conservation(iterations,step=0):
    
    mH = []
    
    while iterations > step :
    
        nx, ny, nz,u,v,H,dt,dx,dy,dz = rand_initialize()
        n = nx * ny
        sum_H0 = 0

        for i in range (1,nx):
            for j in range (1,ny):
                sum_H0 += H[i,j]

        H = fn.H_time_step(H,u,v,nx,ny,dt)  
            
        sum_H1 = 0
    
        for i in range (1,nx):
            for j in range (1,ny):
                sum_H1 += H[i,j]
    

        mH.append(abs((sum_H0-sum_H1)/n))
    
        step+=1

    d_eta_max = max(mH) 
    assert 0.0000000001 > d_eta_max, "mass conservation failed"
    print("mass conservation verifed")


#to verify the momentum conservation in the advection parallel terms:


def momentum_conservation(iterations,step=0):
    
    mv = []
    mu = []    
    
    while iterations > step :
        


        nx, ny, nz,u,v,H,dt,dx,dy,dz = rand_initialize()
        n = nx * ny

        udxu = fn.udexu(u, nx, ny, nz, dx)
        vdyv = fn.vdeyv(v, nx, ny, nz, dy)
        sum_adv_u=0
        sum_adv_v=0
        
        for i in range (1,nx-1):
            for j in range (1,ny-1):
                for k in range (0,2):
                    
                    sum_adv_u+=udxu[k,i,j]
                    sum_adv_v+=vdyv[k,i,j]


        mu.append(sum_adv_u/n)
        mv.append(sum_adv_v/n)

        step+=1
    
    
    d_u_max = max(mu)
    d_v_max = max(mv)

        
    assert 0.0000000001 > d_u_max, "x-momentum conservation failed"
    assert 0.0000000001 > d_v_max, "y-momentum conservation failed"
    
    print("momentum conservation verifed")




def verify_BC_evolution(iterations,step=0):
    while iterations > step :
    
        nx, ny, nz,u,v,H,dt,dx,dy,dz = rand_initialize()
    
        Fx = np.zeros(((nz,nx,ny)))
        Fy = np.zeros(((nz,nx,ny)))
        
        u,v,H,udiff,vdiff,Hdiff = fn.vel_time_step(u, v, H, Fx, Fy, dt, nx, ny)
       
        u_abs1 = abs(u[:,0,:])
        u_abs2 = abs(u[:,-1,:])
        u_abs3 = abs(u[:,:,0])
        u_abs4 = abs(u[:,:,-1])
        
        v_abs1 = abs(v[:,0,:])
        v_abs2 = abs(v[:,-1,:])
        v_abs3 = abs(v[:,:,0])
        v_abs4 = abs(v[:,:,-1])
        
        
        uabs1 = u_abs1.flatten()
        uabs2 = u_abs2.flatten()
        uabs3 = u_abs3.flatten()
        uabs4 = u_abs4.flatten()
        
        vabs1 = v_abs1.flatten()
        vabs2 = v_abs2.flatten()
        vabs3 = v_abs3.flatten()
        vabs4 = v_abs4.flatten()
        
        
        
        Hdiff1 = abs(H[0,1:-1] - H[1,1:-1])
        Hdiff2 = abs(H[-1,1:-1]- H[-2,1:-1])
        Hdiff3 = abs(H[1:-1,0] - H[1:-1,1])
        Hdiff4 = abs(H[1:-1,-1] - H[1:-1,-2])
        
        
        
        umax1 = max(uabs1)
        umax2 = max(uabs2)
        umax3 = max(uabs3)
        umax4 = max(uabs4)
        
        vmax1 = max(vabs1)
        vmax2 = max(vabs2)
        vmax3 = max(vabs3)
        vmax4 = max(vabs4)
        
        Hdmax1 = max(Hdiff1)
        Hdmax2 = max(Hdiff2)
        Hdmax3 = max(Hdiff3)
        Hdmax4 = max(Hdiff4)
        
        assert Hdmax1 == 0, "H BC evolution failed"
        assert Hdmax2 == 0, "H BC evolution failed"
        assert Hdmax3 == 0, "H BC evolution failed"
        assert Hdmax4 == 0, "H BC evolution failed"

        assert umax1 == 0, "u BC evolution failed"
        assert umax2 == 0, "u BC evolution failed"
        assert umax3 == 0, "u BC evolution failed"
        assert umax4 == 0, "u BC evolution failed"

        assert vmax1 == 0, "v BC evolution failed"
        assert vmax2 == 0, "v BC evolution failed"
        assert vmax3 == 0, "v BC evolution failed"
        assert vmax4 == 0, "v BC evolution failed"

        step += 1
    
    print("BC evolution verified")

    
                
        


