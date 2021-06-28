import numpy
from configparser import ConfigParser
import testing_functions as test

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


datafile_u_top="initial_conditions/u_top.txt"
datafile_u_bot="initial_conditions/u_bot.txt"
datafile_v_top="initial_conditions/v_top.txt"
datafile_v_bot="initial_conditions/v_bot.txt"
datafile_u_wind="initial_conditions/u_wind.txt"
datafile_v_wind="initial_conditions/v_wind.txt"
datafile_eta= "initial_conditions/eta.txt"

data_u_top = numpy.loadtxt(datafile_u_top).T
data_v_top = numpy.loadtxt(datafile_v_top).T
data_u_bot = numpy.loadtxt(datafile_u_bot).T
data_v_bot = numpy.loadtxt(datafile_v_bot).T
H = numpy.loadtxt(datafile_eta).T
uw = numpy.loadtxt(datafile_u_wind).T
vw = numpy.loadtxt(datafile_v_wind).T

nx_u_top =len(data_u_top[0])
ny_u_top =len(data_u_top[:,0])

nx_v_top =len(data_v_top[0])
ny_v_top =len(data_v_top[:,0])

nx_u_bot =len(data_u_bot[0])
ny_u_bot =len(data_u_bot[:,0])

nx_v_bot =len(data_v_bot[0])
ny_v_bot =len(data_v_bot[:,0])

nx_eta = len(H[0])
ny_eta = len(H[:,0])

u = numpy.zeros(((2,nx,ny)))

u[1,:,:] = data_u_top[:,:]
u[0,:,:] = data_u_bot[:,:]
a = len(u[:,0,0])
print(a)


v = numpy.zeros(((2,nx,ny)))

v[1,:,:] = data_v_top[:,:]
v[0,:,:] = data_v_bot[:,:]

#create grid



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

test.input(nx,ny,H,nx_u_top,nx_v_top,nx_u_bot,nx_v_bot,nx_eta,ny_u_top,ny_v_top,ny_u_bot,ny_v_bot,ny_eta,uw,vw,u,v,x,y,z,fco,g,nu)
