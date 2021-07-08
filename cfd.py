import numpy
import functions as fn
import input_testing_functions as test
from configparser import ConfigParser

#read const.txt

parser = ConfigParser()
parser.read('initial_conditions/const.txt')
nx = parser.getint('numerical_variables', 'nx') 
ny = parser.getint('numerical_variables', 'ny') 
x = parser.getfloat('physical_variables', 'x') 
y = parser.getfloat('physical_variables', 'y') 
z = parser.getfloat('physical_variables', 'z') 
fco = parser.getfloat('physical_variables', 'fco',
                         fallback = 0) 
g = parser.getfloat('physical_variables', 'g',
                         fallback = 9.81) 
nu = parser.getfloat('physical_variables', 'nu',
                         fallback = 0.004) 

#read IC



datafile_u_wind="initial_conditions/u_wind.txt"
datafile_v_wind="initial_conditions/v_wind.txt"


uw = numpy.loadtxt(datafile_u_wind).T
vw = numpy.loadtxt(datafile_v_wind).T

nx_u_wind =len(uw[0])
ny_u_wind =len(uw[:,0])

nx_v_wind =len(vw[0])
ny_v_wind =len(vw[:,0])

u = numpy.zeros(((2,nx,ny)))
v = numpy.zeros(((2,nx,ny)))
H =  numpy.zeros((nx+1,ny+1))

#create grid

dt = 0.9
nz = 2
dx = x / (nx - 1)
dy = y / (ny - 1)
dz = z / (nz - 1)
xc = numpy.linspace(0, x, nx)
yc = numpy.linspace(0, y, ny)
zc = numpy.linspace(0, z, nz)
X, Y = numpy.meshgrid(xc, yc)

xp = numpy.linspace(0, x, nx+1)
yp = numpy.linspace(0, y, ny+1)
zp = numpy.linspace(0, z, nz)
Xp, Yp = numpy.meshgrid(xp, yp)

#calculate wind stress and bottom stress

Fx,Fy = fn.wind_stress(uw, vw)

#find stationary solution
u, v, H, udiff0, vdiff0, Hdiff0  = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt ,g,fco,nu)

diff0 = udiff0 + vdiff0 + Hdiff0
a = diff0/10000
b = diff0/10000
c = diff0/10000

udiff = 1.1
vdiff = 1.1
Hdiff = 1.1
stepcount = 0
while  udiff > a or vdiff > b or Hdiff > c :

    u, v, H, udiff, vdiff, Hdiff  = fn.vel_time_step(u,v,z,H,Fx,Fy,dx,dy,dz ,dt ,g,fco,nu)
    test.test_eta_H(z,H, nx, ny)
    stepcount +=1

utop= numpy.zeros((nx,ny))
vtop= numpy.zeros((nx,ny))
ubot= numpy.zeros((nx,ny))
vbot= numpy.zeros((nx,ny))

utop[:,:]=u[1,:,:]
vtop[:,:]=v[1,:,:]
ubot[:,:]=u[0,:,:]
vbot[:,:]=v[0,:,:]

#write on file

numpy.savetxt('output_data/uwind', uw)
numpy.savetxt('output_data/vwind', vw)
numpy.savetxt('output_data/utop_out', utop)
numpy.savetxt('output_data/vtop_out', vtop)
numpy.savetxt('output_data/ubot_out', ubot)
numpy.savetxt('output_data/vbot_out', vbot)
numpy.savetxt('output_data/eta_out', H)




