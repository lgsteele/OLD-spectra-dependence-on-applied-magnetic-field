import numpy as np
from numpy import linalg as LA # used in eigenvalue function
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import datetime

##from NVeigenvalues import eigenvalues
##from lor8 import lor8
from DEwidthZeroFieldFit import de


# Generate phase space array of theta, phi, and Bmag
def s2c(theta,phi,Bmag):
    thetaphi = np.array([[x0,y0] for x0 in theta \
                      for y0 in phi])
    theta,phi = thetaphi.transpose()
    xyz = (phi*np.array([np.cos(phi)*np.sin(theta),\
               np.sin(phi)*np.sin(theta),\
               np.cos(theta)])).transpose()
    Bxyztmp = (Bmag*xyz).reshape(len(Bmag),len(thetaphi),1,3)
    return Bxyztmp


# Calculate the eigenvalues for phase space
def eigenvalues(Bxyz):
    Bxyztmp = np.copy(Bxyz)
    base = np.ones((len(Bxyztmp)*len(Bxyztmp[0]),1)).\
           reshape(len(Bxyztmp),len(Bxyztmp[0]),1,1)
    zf = np.array([2.87e9,3.6e6,2.6e6,0,0,0,0,0])
    zfs1 = zf[0]*(base*np.array([[0,.3333,.3333],\
                    [.3333,0,-.3333],\
                    [.3333,-.3333,0]]))
    strain1 = zf[1]*(base*np.array([[.3333,-.3333,.6667],\
                    [-.3333,-.6667,.3333],\
                    [.6667,.3333,.3333]]))
    zfs2 = zf[0]*(base*np.array([[0,.3333j,-.3333],\
                    [-.3333j,0,-.3333j],\
                    [-.3333,.3333j,0]]))
    strain2 = zf[1]*(base*np.array([[.3333,-.3333j,-.6667],\
                    [.3333j,-.6667,.3333j],\
                    [-.6667,-.3333j,.3333]]))
    zfs3 = zf[0]*(base*np.array([[0,-.3333,.3333],\
                    [-.3333,0,.3333],\
                    [.3333,.3333,0]]))
    strain3 = zf[1]*(base*np.array([[.3333,.3333,.6667],\
                    [.3333,-.6667,-.3333],\
                    [.6667,-.3333,.3333]]))
    zfs4 = zf[0]*(base*np.array([[0,-.3333j,-.3333],\
                    [.3333j,0,.3333j],\
                    [-.3333,-.3333j,0]]))
    strain4 = zf[1]*(base*np.array([[.3333,.3333j,-.6667],\
                    [-.3333j,-.6667,-.3333j],\
                    [-.6667,.3333j,.3333]]))
    sx = base*np.array([[0,.7071,0],[.7071,0,.7071],[0,.7071,0]])
    sy = base*np.array([[0,-.7071j,0],[-.7071j,0,-.7071j],[0,-.7071j,0]])
    sz = base*np.array([[1,0,0],[0,0,0],[0,0,-1]])
    zeeman = 28024951642 * ((Bxyztmp[:,:,:,0].reshape(len(Bxyztmp),\
                                        len(Bxyztmp[0]),1,1)*sx) + \
                            (Bxyztmp[:,:,:,1].reshape(len(Bxyztmp),\
                                        len(Bxyztmp[0]),1,1)*sy) + \
                            (Bxyztmp[:,:,:,2].reshape(len(Bxyztmp),\
                                        len(Bxyztmp[0]),1,1)*sz))
    w, v = LA.eigh(zfs1+strain1+zeeman)
    f1 = (abs(w[:,:,0])+abs(w[:,:,1]))
    f2 = (abs(w[:,:,0])+abs(w[:,:,2]))
    w, v = LA.eigh(zfs2+strain2+zeeman)
    f3 = (abs(w[:,:,0])+abs(w[:,:,1]))
    f4 = (abs(w[:,:,0])+abs(w[:,:,2]))
    w, v = LA.eigh(zfs3+strain3+zeeman)
    f5 = (abs(w[:,:,0])+abs(w[:,:,1]))
    f6 = (abs(w[:,:,0])+abs(w[:,:,2]))
    w, v = LA.eigh(zfs4+strain4+zeeman)
    f7 = (abs(w[:,:,0])+abs(w[:,:,1]))
    f8 = (abs(w[:,:,0])+abs(w[:,:,2]))
    evals = np.dstack([f1,f2,f3,f4,f5,f6,f7,f8])
    evals = np.sort(evals,axis=2).reshape(len(Bxyztmp),\
                                len(Bxyztmp[0]),1,8)
##    return np.concatenate((Bxyztmp,evals),axis=3)
    return evals


# Generate spectra in phase space
def lor8(freq,zf,ev):
    evtmp = np.copy(ev)
    base = np.ones((len(evtmp)*len(evtmp[0]),1)).\
           reshape(len(evtmp),len(evtmp[0]),1,1)
    freq = freq*base
    ev1 = (ev[:,:,:,0].reshape(len(evtmp),len(evtmp[0]),1,1))\
            *base
    ev2 = (ev[:,:,:,1].reshape(len(evtmp),len(evtmp[0]),1,1))\
            *base
    ev3 = (ev[:,:,:,2].reshape(len(evtmp),len(evtmp[0]),1,1))\
            *base
    ev4 = (ev[:,:,:,3].reshape(len(evtmp),len(evtmp[0]),1,1))\
            *base
    ev5 = (ev[:,:,:,4].reshape(len(evtmp),len(evtmp[0]),1,1))\
            *base
    ev6 = (ev[:,:,:,5].reshape(len(evtmp),len(evtmp[0]),1,1))\
            *base
    ev7 = (ev[:,:,:,6].reshape(len(evtmp),len(evtmp[0]),1,1))\
            *base
    ev8 = (ev[:,:,:,7].reshape(len(evtmp),len(evtmp[0]),1,1))\
            *base
    lor1 = 1e7*((zf[2]**2)/\
            (np.pi*zf[2]*((freq-ev1)**2+ zf[2]**2)))
    lor2 = 1e7*((zf[2]**2)/\
            (np.pi*zf[2]*((freq-ev2)**2+ zf[2]**2)))
    lor3 = 1e7*((zf[2]**2)/\
            (np.pi*zf[2]*((freq-ev3)**2+ zf[2]**2)))
    lor4 = 1e7*((zf[2]**2)/\
            (np.pi*zf[2]*((freq-ev4)**2+ zf[2]**2)))
    lor5 = 1e7*((zf[2]**2)/\
            (np.pi*zf[2]*((freq-ev5)**2+ zf[2]**2)))
    lor6 = 1e7*((zf[2]**2)/\
            (np.pi*zf[2]*((freq-ev6)**2+ zf[2]**2)))
    lor7 = 1e7*((zf[2]**2)/\
            (np.pi*zf[2]*((freq-ev7)**2+ zf[2]**2)))
    lor8 = 1e7*((zf[2]**2)/\
            (np.pi*zf[2]*((freq-ev8)**2+ zf[2]**2)))
    return lor1+lor2+lor3+lor4+lor5+lor6+lor7+lor8


# Fit spectra and return splitting and widths



############################################
# Define theta, phi, Bmag arrays to calculate eigenvalues
############################################
##theta = np.linspace(0.01,np.pi/2.,num=2)
##phi = np.linspace(0.,np.pi/2.,num=2)
theta = np.array([1,2])
phi = np.array([1,2])
Bmag = np.linspace(1e-7,1e-6,2)
Bmag = Bmag.reshape(len(Bmag),1,1)
Bxyz = s2c(theta,phi,Bmag)
ev = eigenvalues(Bxyz)

############################################
# Define frequency range and zf for generating spectra
############################################
freq = np.linspace(2.77e9,2.97e9,1e6)
zf = np.array([2.87e9,3.6e6,2.6e6,0,0,0,0,0])
spectra = lor8(freq,zf,ev)
print freq.shape
print spectra.shape
plt.plot(freq,spectra[0,0,0,:])
plt.show()
##NVspectrafromB(freq,ev)


############################################
############################################








