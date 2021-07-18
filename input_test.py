import numpy as np
from configparser import ConfigParser

#tests to do list

def test_grid(nx,ny,nx_uw,nx_vw,ny_uw,ny_vw):
    v_nx_uw = nx-nx_uw
    v_nx_vw = nx-nx_vw
    v_ny_uw = ny-ny_uw
    v_ny_vw = ny-ny_vw

    assert nx > 2
    
    assert ny > 2
    
    assert v_nx_uw == 0, "x dimension of u_top doesn't match with the declared grid"
    
    assert v_nx_vw == 0, "x dimension of v_top doesn't match with the declared grid" 
    
    assert v_ny_uw == 0, "y dimension of u_top doesn't match with the declared grid"
    
    assert v_ny_vw == 0, "y dimension of v_top doesn't match with the declared grid" 
     


def test_input_values(nx,ny,uw,vw):
    
    for i in range (0,nx):
        for j in range (0,ny):
                
            a_uw = abs(uw[i,j])
            a_vw = abs(vw[i,j])
                
            assert a_uw < 30, "please provide zonal wind values below 30 m/s"
            assert a_vw < 30, "please provide meridional wind values below 30 m/s"
                
            
def test_ph_parameters(x,y,z,fco,g,nu):
    afco = abs(fco)
    assert x > 0, "please provide positive dimensional values for the basin"
    assert y > 0, "please provide positive dimensional values for the basin"
    assert z > 0, "please provide positive dimensional values for the basin"
    assert 0.2 > afco, "please provide a physically plausible value for coriolis parameter (around ±0.0001 values)"
    assert g > 0, "gravitational acceleration [g] must be positive! (positive values are associated to downward vectors)"
    assert nu > 0, "viscosity [nu] must be positive!"

def input_test(nx,ny,uw,vw,nx_uw,nx_vw,ny_uw,ny_vw,x,y,z,fco,g,nu)  :
    
    test_input_values(nx,ny,uw,vw)
    test_grid(nx,ny,nx_uw,nx_vw,ny_uw,ny_vw)
    test_ph_parameters(x,y,z,fco,g,nu)
    print("input values verified")  

#read const.txt

parser = ConfigParser()
parser.read('initial_conditions/const.txt')
nx = parser.getint('numerical_variables', 'nx',
                         fallback = -1) 
ny = parser.getint('numerical_variables', 'ny',
                         fallback = -1) 

nz = 2
x = parser.getfloat('physical_variables', 'x',
                         fallback = -1) 
y = parser.getfloat('physical_variables', 'y',
                         fallback = -1) 
z = parser.getfloat('physical_variables', 'z',
                         fallback = -1) 
fco = parser.getfloat('physical_variables', 'fco',
                         fallback = -1) 
g = parser.getfloat('physical_variables', 'g',
                         fallback = -1) 
nu = parser.getfloat('physical_variables', 'nu',
                         fallback = -1) 

#read IC



datafile_u_wind="initial_conditions/u_wind.txt"
datafile_v_wind="initial_conditions/v_wind.txt"




uw = np.loadtxt(datafile_u_wind).T
vw = np.loadtxt(datafile_v_wind).T

ny_uw =len(uw[0])
nx_uw =len(uw[:,0])

ny_vw =len(vw[0])
nx_vw =len(vw[:,0])

#create grid

dx = x / (nx - 1)
dy = y / (ny - 1)
dz = z / (nz - 1)
xc = np.linspace(0, x, nx)
yc = np.linspace(0, y, ny)
zc = np.linspace(0, z, nz)
X, Y = np.meshgrid(xc, yc)

xp = np.linspace(0, x, nx+1)
yp = np.linspace(0, y, ny+1)
zp = np.linspace(0, z, nz)
Xp, Yp = np.meshgrid(xp, yp)

#apply the tests

test_input_values(nx,ny,uw,vw)
print("input values verified")
test_grid(nx,ny,nx_uw,nx_vw,ny_uw,ny_vw)
print("grid dimensions verified")
test_ph_parameters(x,y,z,fco,g,nu)
print("physical parameters verified")

print("######################")
print("#  ALL TESTS PASSED  #")
print("######################")
