import numpy as np
import functions as fn

#input test

def test_grid(nx,ny,nx_u_top,nx_v_top,nx_u_bot,nx_v_bot,nx_eta,ny_u_top,ny_v_top,ny_u_bot,ny_v_bot,ny_eta):
    v_nx_u_top = nx-nx_u_top
    v_nx_v_top = nx-nx_v_top
    v_nx_u_bot = nx-nx_u_bot
    v_nx_v_bot = nx-nx_v_bot
    v_nx_eta = nx+1-nx_eta
    v_ny_u_top = ny-ny_u_top
    v_ny_v_top = ny-ny_v_top
    v_ny_u_bot = ny-ny_u_bot
    v_ny_v_bot = ny-ny_v_bot
    v_ny_eta = ny+1-ny_eta
    
    
    assert nx > 2
    
    assert ny > 2
    
    assert v_nx_u_top == 0, "x dimension of u_top doesn't match with the declared grid"
    
    assert v_nx_v_top == 0, "x dimension of v_top doesn't match with the declared grid" 
    
    assert v_nx_u_bot == 0, "x dimension of u_bot doesn't match with the declared grid"
    
    assert v_nx_v_bot == 0, "x dimension of v_bot.txt doesn't match with the declared grid"
    
    assert v_nx_eta == 0, "x dimension of eta.txt doesn't match with the declared grid"
    
    assert v_ny_u_top == 0, "y dimension of u_top doesn't match with the declared grid"
    
    assert v_ny_v_top == 0, "y dimension of v_top doesn't match with the declared grid" 
    
    assert v_ny_u_bot == 0, "y dimension of u_bot doesn't match with the declared grid"
    
    assert v_ny_v_bot == 0, "y dimension of v_bot.txt doesn't match with the declared grid"
    
    assert v_ny_eta == 0, "y dimension of eta.txt doesn't match with the declared grid"
    


def test_input_values(nx,ny,u,v,uw,vw):
    
    for i in range (0,nx):
        for j in range (0,ny):
            for k in range (0,2):
                a_u = abs(u[k,i,j])
                a_v = abs(v[k,i,j])
                a_uw = abs(u[k,i,j])
                a_vw = abs(v[k,i,j])
                
                assert a_u < 1.1, "please provide zonal current values below 1.1 m/s"
                assert a_v < 1.1, "please provide meridional current values below 1.1 m/s"
                assert a_uw < 30, "please provide zonal wind values below 30 m/s"
                assert a_vw < 30, "please provide meridional wind values below 30 m/s"
                
            

def test_eta_H(z,H, nx, ny):
    
    for i in range (0,nx):
        for j in range (0,ny):
            a_H = abs(H[i,j])
            z_half = z/2
            diff = z_half-a_H
            
            assert diff  > 0, "ERROR : h/2+eta < 0, non-physical solution, try higher values of z"
                
def test_eta_mean_null(H, nx, ny):
   
    a=0
    for i in range (1,nx):
        for j in range (1,ny):
            a+=H(i,j)
    
    assert a == 0, "please provide basin surface elevation values with mean 0"
    
def test_input_eta_values(nx,ny,H):
    
    for i in range (0,nx+1):
        for j in range (0,ny+1):
            a = abs(H[i,j])
            
            assert a < 2, "please provide basin surface elevation values below 2 m"
            
def test_vel_x_boundaries(ny,u,v):
    
    for j in range (0,ny):
        for k in range (0,2):
            a = u[k,0,j]
            b = u[k,-1,j]
            c = v[k,0,j]
            d = v[k,-1,j]
            
            assert a  == 0, "please provide null zonal velocity values at the boundary x=0"
            assert b  == 0, "please provide null zonal velocity values at the boundary x=L"
            assert c  == 0, "please provide null meridional velocity values at the boundary x=0"
            assert d  == 0, "please provide null meridional velocity values at the boundary x=L"
            
def test_vel_y_boundaries(nx,u,v):
    
    for i in range (0,nx):
        for k in range (0,2):
            a = u[k,i,0]
            b = u[k,i,-1]
            c = v[k,i,0]
            d = v[k,i,-1]
            
            assert a  == 0, "please provide null zonal velocity values at the boundary y=0"
            assert b  == 0, "please provide null zonal velocity values at the boundary y=L"
            assert c  == 0, "please provide null meridional velocity values at the boundary y=0"
            assert d  == 0, "please provide null meridional velocity values at the boundary y=L"
            
def test_eta_x_boundaries(ny,H):
    
    for j in range (1,ny):
        a0 = H[0,j]
        a1 = H[1,j]
        b0 = H[-1,j]
        b1 = H[-2,j]
        a = a0-a1
        b = b0-b1
        assert a  == 0, "please provide null elevation gradient values at the boundary x=0"
        assert b  == 0, "please provide null elevation gradient values at the boundary x=L"
        
        
        
def test_eta_y_boundaries(nx,H):
    
    for i in range (1,nx):
        a0 = H[i,0]
        a1 = H[i,1]
        b0 = H[i,-1]
        b1 = H[i,-2]
        a = a0-a1
        b = b0-b1
        assert a  == 0, "please provide null elevation gradient values at the boundary y=0"
        assert b  == 0, "please provide null elevation gradient values at the boundary y=L"      
 

def test_ph_parameters(x,y,z,fco,g,nu):
    afco = abs(fco)
    assert x > 0, "please provide positive dimensional values for the basin"
    assert y > 0, "please provide positive dimensional values for the basin"
    assert z > 0, "please provide positive dimensional values for the basin"
    assert 0.001 > afco, "please provide a physically plausible value for coriolis parameter (around ±0.0001 values)"
    assert g > 0, "gravitational acceleration [g] must be positive! (positive values are associated to downward vectors)"
    assert nu > 0, "viscosity [nu] must be positive!"

def input(nx,ny,H,nx_u_top,nx_v_top,nx_u_bot,nx_v_bot,nx_eta,ny_u_top,ny_v_top,ny_u_bot,ny_v_bot,ny_eta,uw,vw,u,v,x,y,z,fco,g,nu)  :
    test_eta_y_boundaries(nx,H)
    test_eta_x_boundaries(ny,H)
    test_vel_y_boundaries(nx,u,v)
    test_vel_x_boundaries(ny,u,v)
    test_input_eta_values
    test_eta_mean_null
    test_eta_H    
    test_input_values(nx,ny,u,v,uw,vw)
    test_grid(nx,ny,nx_u_top,nx_v_top,nx_u_bot,nx_v_bot,nx_eta,ny_u_top,ny_v_top,ny_u_bot,ny_v_bot,ny_eta)
    test_ph_parameters(x,y,z,fco,g,nu)
    print("input values verified")  