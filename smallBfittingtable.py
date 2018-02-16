import numpy as np
from numpy import linalg as LA # used in eigenvalue function
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import datetime

from NVeigenvalues import eigenvalues
from lor8 import lor8
from DEwidthZeroFieldFit import de


# Generate Bxyz array from theta, phi, a Bmag
def s2c(theta,phi,Bmag):
    thetaphi = np.array([[x0,y0] for x0 in theta \
                      for y0 in phi])
    theta,phi = thetaphi.transpose()
    xyz = (phi*np.array([np.cos(phi)*np.sin(theta),\
               np.sin(phi)*np.sin(theta),\
               np.cos(theta)])).transpose()
    Bxyztmp = (Bmag*xyz).reshape(2,4,1,3)
    return Bxyztmp


# Calculate the eigenvalues for array of B-field
# magnitudes and orientations as defined by
# theta, phi, and Bmag
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
    return np.concatenate((Bxyztmp,evals),axis=3)

############################################
# Next is to generate spectra for each Bxyz...
def NVspectrafromB(freq,zf,amp,Bx,By,Bz):
    # Define relevant arrays
    # Calculate eigenvalues for given Bx,By,Bz
    ev = eigenvalues(zf[0:3],[Bx,By,Bz])
    # Generate NV- spectra
    spectra = lor8(freq,zf,amp,ev)
##    plt.plot(freq,spectra)
##    plt.show()
    return spectra
##NVspectrafromB()


############################################
# Define theta, phi, Bmag arrays to calculate eigenvalues
############################################
##theta = np.linspace(0.01,np.pi/2.,num=2)
##phi = np.linspace(0.,np.pi/2.,num=2)
theta = np.array([1,2])
phi = np.array([1,2])
Bmag = np.linspace(1e-7,1e-6,2).reshape(2,1,1)
Bxyz = s2c(theta,phi,Bmag)
ev = eigenvalues(Bxyz)
print ev

############################################
# Define frequency range for generating spectra
############################################
freq = np.linspace(2.77e9,2.97e9,1e6)



############################################
############################################








