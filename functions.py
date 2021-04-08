import numpy
from configparser import ConfigParser

# Choice of physical an numerical constants 
parser = ConfigParser()
parser.read('const.txt')

x = parser.getfloat('physical_variables', 'x',
                         fallback = 10)
y = parser.getfloat('physical_variables', 'y',
                         fallback = 10)
z = parser.getfloat('physical_variables', 'z',
                         fallback = 1) 
g = parser.getfloat('physical_variables', 'g',
                         fallback = 10) 
fco = parser.getfloat('physical_variables', 'fco',
                          fallback = 0)
nu = parser.getfloat('physical_variables', 'nu',
                          fallback = 0)
nx = parser.getfloat('numerical_variables', 'nx',
                         fallback = -1) 
ny = parser.getfloat('numerical_variables', 'ny',
                         fallback = -1) 
nz = 2
dx = x / (nx - 1)
dy = y / (ny - 1)
dz = z / (nz - 1)

#calcola stress del fondale

def bottom_stress(u, v, nx, ny, nz):
    Bx = numpy.zeros(((nz,nx,ny)))
    By = numpy.zeros(((nz,nx,ny)))
    k = 0.01
    Bx[0,:,:]= -k*u[0,:,:]*numpy.sqrt((u[0,:,:]**2)+(v[0,:,:]**2))
    By[0,:,:]= -k*v[0,:,:]*numpy.sqrt((u[0,:,:]**2)+(v[0,:,:]**2))
    return Bx, By

#calcola stress ventoso

def wind_stress(uw, vw, nx, ny, nz):
    Fx = numpy.zeros(((nz,nx,ny)))
    Fy = numpy.zeros(((nz,nx,ny)))
    k = 0.001
    Fx[1,:,:]= k*uw[:,:]*numpy.sqrt((uw[:,:]**2)+(vw[:,:]**2))
    Fy[1,:,:]= k*vw[:,:]*numpy.sqrt((uw[:,:]**2)+(vw[:,:]**2))
    return Fx, Fy

#derivata in x funzione in due variabili

def Dexb(f, nx, ny, dx):
    f_1= numpy.zeros((nx+1,ny+1))
    f_1[:-1,:]=(f[1:,:]-f[:-1,:])/dx
    return f_1

#derivata in y funzione in due variabili

def Deyb(f, nx, ny, dy):
    f_1 =  numpy.zeros((nx+1,ny+1))
    f_1[:,:-1] = (f[:,1:]-f[:,:-1])/dy
    return f_1

#derivata in x funzione in tre variabili

def Dex(f, nx, ny, nz, dx):
    f_1 = numpy.zeros(((nz,nx,ny)))
    f_1[:,:-1,:] = (f[:,1:,:]-f[:,:-1,:])/dx
    return f_1

#derivata in y funzione in tre variabili

def Dey(f, nx, ny, nz, dy):
    f_1 = numpy.zeros(((nz,nx,ny)))
    f_1[:,:,:-1] = (f[:,:,1:]-f[:,:,:-1])/dy
    return f_1

#derivata seconda in x funzione in tre variabili

def Dex2(f, nx, ny, nz, dx):
    f_2 = numpy.zeros(((nz,nx,ny)))
    f_2[:,1:-1,:] = (f[:,2:nx,:]+f[:,0:nx-2,:]-(2*f[:,1:-1,:]))/dx**2
    return f_2

#derivata seconda in y funzione in tre variabili

def Dey2(f, nx, ny, nz, dy):
    f_2 = numpy.zeros(((nz,nx,ny)))
    f_2[:,:,1:-1] = (f[:,:,2:ny]+f[:,:,0:ny-2]-(2*f[:,:,1:-1]))/dy**2
    return f_2

#derivata seconda in x funzione in tre variabili
def Dex2b(f):
    f_2 = numpy.zeros((nx+1,ny+1))
    f_2[1:-1,:] = (f[2:,:]+f[0:-2,:]-(2*f[1:-1,:]))/dx**2
    return f_2

#derivata seconda in y funzione in tre variabili

def Dey2b(f, nx, ny, dy):
    f_2 = numpy.zeros((nx+1,ny+1))
    f_2[:,1:-1] = (f[:,2:]+f[:,0:-2]-(2*f[:,1:-1]))/dy**2
    return f_2

def H_time_step(H,u,v,nx,ny):
    Hn = H.copy()
    U= numpy.zeros((nx+1,ny+1))
    V= numpy.zeros((nx+1,ny+1))
    U[1:,1:]=(sum(u[:,:,:]))*(z+Hn[:-1,:-1])/nz
    V[1:,1:]=(sum(v[:,:,:]))*(z+Hn[:-1,:-1])/nz
   
    DexbU = Dexb(U,nx,ny,dx)
    DeybV = Deyb(V,nx,ny,dy)
    H[1:-1,1:-1]=Hn[1:-1,1:-1]-(DexbU[1:-1,1:-1]+DexbU[1:-1,2:])/2-(DeybV[1:-1,1:-1]+DeybV[2:,1:-1])/2
    #BC gradiente di pressione nullo al bordo lungo la perpendicolare
    H[:,0] = H[:,1]
    H[:,ny]=H[:,ny-1]
    H[0,:]  = H[1,:]
    H[nx,:] = H[nx-1,:]

    return H


def vel_time_step(u,v,H,Fx,Fy,dt, nx, ny):
    Hn = H.copy()
    H = H_time_step(H,u,v,nx,ny)
    
    Bx,By = bottom_stress(u, v, nx, ny, nz)
    
    cox = numpy.zeros(((nz,nx,ny)))
    coy = numpy.zeros(((nz,nx,ny)))
    udxu = numpy.zeros(((nz,nx,ny)))
    udxv = numpy.zeros(((nz,nx,ny)))
    vdyu = numpy.zeros(((nz,nx,ny)))
    vdyv = numpy.zeros(((nz,nx,ny)))
    dexP = numpy.zeros((nx,ny))
    deyP = numpy.zeros((nx,ny))

    disu = numpy.zeros(((nz,nx,ny)))
    disv = numpy.zeros(((nz,nx,ny)))
    Dez2un = numpy.zeros(((nz,nx,ny)))
    Dez2vn = numpy.zeros(((nz,nx,ny)))
    
    un = u.copy()
    vn = v.copy()
    Dexun = Dex(un,nx,ny,nz,dx)
    Dexvn = Dex(vn,nx,ny,nz,dx)
    Deyun = Dey(un,nx,ny,nz,dy)
    Deyvn = Dey(vn,nx,ny,nz,dy)
    Dez2un[0,:,:]=-(un[0,:,:]-un[1,:,:])/(dz**2)
    Dez2un[1,:,:]=-Dez2un[0,:,:]
    Dez2vn[0,:,:]=-(vn[0,:,:]-vn[1,:,:])/(dz**2)
    Dez2vn[1,:,:]=-Dez2vn[0,:,:]
     
    
    cox[:,:,:] = fco*vn[:,:,:]
    coy[:,:,:] = -fco*un[:,:,:]
    udxu[:,1:-1,1:-1] = un[:,1:-1,1:-1]*(Dexun[:,0:-2,1:-1]+Dexun[:,1:-1,1:-1])/2
    udxv[:,1:-1,1:-1] = un[:,1:-1,1:-1]*(Dexvn[:,0:-2,1:-1]+Dexvn[:,1:-1,1:-1])/2
    vdyu[:,1:-1,1:-1] = vn[:,1:-1,1:-1]*(Deyun[:,1:-1,0:-2]+Deyun[:,1:-1,1:-1])/2
    vdyv[:,1:-1,1:-1] = vn[:,1:-1,1:-1]*(Deyvn[:,1:-1,0:-2]+Deyvn[:,1:-1,1:-1])/2
    dexP[:,:] = g/2 * (Dexb(H,nx,ny,dx)[:-1,:-1]+Dexb(H,nx,ny,dx)[:-1,1:])
    deyP[:,:] = g/2 * (Deyb(H,nx,ny,dy)[:-1,:-1]+Deyb(H,nx,ny,dy)[1:,:-1])
    disuh = nu * (Dex2(un,nx,ny,nz,dx) + Dey2(un,nx,ny,nz,dy))
    disvh = nu * (Dex2(vn,nx,ny,nz,dx) + Dey2(vn,nx,ny,nz,dy))
    disu[:,:,:] = disuh[:,:,:] + Dez2un[:,:,:]
    disv[:,:,:] = disvh[:,:,:] + Dez2vn[:,:,:]
    
    u[:,1:-1,1:-1] = (un[:,1:-1,1:-1] - dexP[1:-1,1:-1]-udxu[:,1:-1,1:-1]-vdyu[:,1:-1,1:-1]+disu[:,1:-1,1:-1]+cox[:,1:-1,1:-1]+Fx[:,1:-1,1:-1]+Bx[:,1:-1,1:-1])*dt
    v[:,1:-1,1:-1] = (vn[:,1:-1,1:-1] - deyP[1:-1,1:-1]-udxv[:,1:-1,1:-1]-vdyv[:,1:-1,1:-1]+disv[:,1:-1,1:-1]+coy[:,1:-1,1:-1]+Fy[:,1:-1,1:-1]+By[:,1:-1,1:-1])*dt
    u[:,0,:] = 0
    u[:,-1,:] = 0
    v[:,0,:] = 0
    v[:,-1,:] = 0
    u[:,:,0] = 0
    u[:,:,-1] = 0
    v[:,:,0] = 0
    v[:,:,-1] = 0
    
    du2 = (u-un)**2
    dv2 = (v-vn)**2
    dH2 = (H-Hn)**2
   
    u2 = u**2
    v2 = v**2
    H2 = H**2

    udiff = numpy.sum(du2)/(numpy.sum(u2)+.0000000001)
    vdiff = numpy.sum(dv2)/(numpy.sum(v2)+.0000000001)
    Hdiff = numpy.sum(dH2)/(numpy.sum(H2)+.0000000001)
    
    return u,v,H,udiff,vdiff,Hdiff


