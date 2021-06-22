import numpy
import functions as fn
import testing_functions as test
from configparser import ConfigParser
from matplotlib import pyplot, cm

#read const.txt

parser = ConfigParser()
parser.read('initial_conditions/const.txt')
nx = parser.getint('numerical_variables', 'nx',
                         fallback = -1) 
ny = parser.getint('numerical_variables', 'ny',
                         fallback = -1) 
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
H =  numpy.loadtxt(datafile_eta).T
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


v = numpy.zeros(((2,nx,ny)))

v[1,:,:] = data_v_top[:,:]
v[0,:,:] = data_v_bot[:,:]



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

Fx,Fy = fn.wind_stress(uw, vw, nx, ny, nz)

#find stationary solution
u, v, H, udiff0, vdiff0, Hdiff0  = fn.vel_time_step(u, v, H, Fx, Fy, dt, nx, ny,g)

diff0 = udiff0 + vdiff0 + Hdiff0
a = diff0/100000
b = diff0/100000
c = diff0/100000

udiff = 1.1
vdiff = 1.1
Hdiff = 1.1
stepcount = 0
while  udiff > a or vdiff > b or Hdiff > c :

    u, v, H, udiff, vdiff, Hdiff  = fn.vel_time_step(u, v, H, Fx, Fy, dt, nx, ny, g)
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

#plot

fig = pyplot.figure(figsize = (11,7), dpi=100)
pyplot.quiver(Y, X, utop, vtop);
fig.savefig('output_figures/vel_top_plot.png')

fig = pyplot.figure(figsize = (11,7), dpi=100)
pyplot.quiver(Y, X, ubot, vbot);
fig.savefig('output_figures/vel_bot_plot.png')

fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')                      
surf = ax.plot_surface(Xp, Yp, H[:], cmap=cm.viridis) 
fig.savefig('output_figures/eta_plot.png')  

print(stepcount,"step")



