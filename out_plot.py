import numpy
from configparser import ConfigParser
from matplotlib import pyplot, cm

datafile_u_wind="output_data/uwind"
datafile_v_wind="output_data/vwind"
datafile_u_top="output_data/utop_out"
datafile_u_bot="output_data/ubot_out"
datafile_v_top="output_data/vtop_out"
datafile_v_bot="output_data/vbot_out"
datafile_eta= "output_data/eta_out"

uw = numpy.loadtxt(datafile_u_wind).T
vw = numpy.loadtxt(datafile_v_wind).T
utop = numpy.loadtxt(datafile_u_top).T
vtop = numpy.loadtxt(datafile_v_top).T
ubot = numpy.loadtxt(datafile_u_bot).T
vbot = numpy.loadtxt(datafile_v_bot).T
H =  numpy.loadtxt(datafile_eta).T

nx = len(utop[0])
ny = len(vtop[:,0])
nz = 2

parser = ConfigParser()
parser.read('initial_conditions/const.txt')
x = parser.getfloat('physical_variables', 'x') 
y = parser.getfloat('physical_variables', 'y') 
z = parser.getfloat('physical_variables', 'z')

xc = numpy.linspace(0, x, nx)
yc = numpy.linspace(0, y, ny)
zc = numpy.linspace(0, z, nz)
X, Y = numpy.meshgrid(xc, yc)

xp = numpy.linspace(0, x, nx+1)
yp = numpy.linspace(0, y, ny+1)
zp = numpy.linspace(0, z, nz)
Xp, Yp = numpy.meshgrid(xp, yp)



#plot

fig = pyplot.figure(figsize = (11,11), dpi=100)
pyplot.quiver(X, Y, uw, vw);
pyplot.title("Wind Velocity Field")
fig.savefig('output_figures/wind_plot.png')

fig = pyplot.figure(figsize = (11,11), dpi=100)
pyplot.quiver(X, Y, utop, vtop);
pyplot.title("Surface Current Velocity Field")
fig.savefig('output_figures/vel_top_plot.png')

fig = pyplot.figure(figsize = (11,11), dpi=100)
pyplot.quiver(X, Y, ubot, vbot);
pyplot.title("Bottom Current Velocity Field")
fig.savefig('output_figures/vel_bot_plot.png')

fig = pyplot.figure(figsize=(11, 11), dpi=100)
ax = fig.gca(projection='3d')                      
surf = ax.plot_surface(Xp, Yp, H[:], cmap=cm.viridis)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('H [m]')
ax.set_title('Surface Elevation')
fig.savefig('output_figures/eta_plot.png')  

