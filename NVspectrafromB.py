import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import datetime

from NVeigenvalues import eigenvalues
from lor8 import lor8
from DEwidthZeroFieldFit import de

# This code solves for the splitting a width of a two-lorentzian fit
# of NV-spectra for a few orientations of the B-field
# For a more robust account on how the splitting and width depend
# on field-orientation, see smallBfittingtable


# Define variable arrays
zf = np.array([2.87e9,3.6e6,2.6e6,0,\
                0,0,0,0])
amp = np.array([1e7,1e7,1e7,1e7,1e7,1e7,1e7,1e7])
freq = np.arange(2.845e9,2.895e9,1e5)

# Given a magnetic field [Bx,By,Bz], generate
# NV- spectra
def NVspectrafromB(freq,zf,amp,Bx,By,Bz):
    # Define relevant arrays
    # Calculate eigenvalues for given Bx,By,Bz
    ev = eigenvalues(zf[0:3],[Bx,By,Bz])
    # Generate NV- spectra
    spectra = lor8(freq,zf,amp,ev)
##    plt.plot(freq,spectra)
##    plt.show()
    return spectra






# Next, generate spectra for various values of B
BArray = np.arange(1e-7,2.5e-4,10e-7)
spectrax = np.zeros(( len(BArray),len(freq)+1 ))
spectrax[:,0] = np.copy(BArray)
spectraz = np.copy(spectrax)
spectraxy = np.copy(spectrax)
spectraxz = np.copy(spectrax)
spectraxyz = np.copy(spectrax)
ewx = np.zeros((len(spectrax),8))
ewz = np.zeros((len(spectrax),8))
##ewxy = np.zeros((len(spectrax),8))
ewxz = np.zeros((len(spectrax),8))
ewxyz = np.zeros((len(spectrax),8))

for i in range(0,len(spectrax),1):
    spectrax[i,1:] = NVspectrafromB(freq,zf,amp,BArray[i],0,0)
    spectraz[i,1:] = NVspectrafromB(freq,zf,amp,0,0,BArray[i])
##    spectraxy[i,1:] = NVspectrafromB(freq,zf,amp,\
##                                     BArray[i]/np.sqrt(2),BArray[i]/np.sqrt(2),0)
    spectraxz[i,1:] = NVspectrafromB(freq,zf,amp,\
                                     BArray[i]/np.sqrt(2),0,BArray[i]/np.sqrt(2))
    spectraxyz[i,1:] = NVspectrafromB(freq,zf,amp,\
                                     BArray[i]/np.sqrt(3),\
                                      BArray[i]/np.sqrt(3),\
                                      BArray[i]/np.sqrt(3))
    ewx[i,:] = np.copy(de(freq,spectrax[i,1:],zf))
    ewz[i,:] = np.copy(de(freq,spectraz[i,1:],zf))
##    ewxy[i,:] = np.copy(de(freq,spectraxy[i,1:],zf))
    ewxz[i,:] = np.copy(de(freq,spectraxz[i,1:],zf))
    ewxyz[i,:] = np.copy(de(freq,spectraxyz[i,1:],zf))
##print ewxy[:,1]
##print len(spectrax)
##print len(spectraz)
##print len(spectraxy)
##print de(freq,spectraxy[i,1:],zf)
##print de(freq,spectraxy[249,1:],zf)




fig = plt.figure(figsize=plt.figaspect(1.))

ax = fig.add_subplot(2,3,1)
X = freq/1e9
Y = spectrax[:,0]*1e6
X,Y = np.meshgrid(X,Y)
Z = spectrax[:,1:]
plt.contourf(X,Y,Z,100,cmap='Reds')
plt.colorbar()
ax.set_title('Bx')
ax.set_xlabel('Freq (GHz)')
ax.set_ylabel('B (uT)')

ax = fig.add_subplot(2,3,2)
X = freq/1e9
Y = spectraz[:,0]*1e6
X,Y = np.meshgrid(X,Y)
Z = spectraz[:,1:]
plt.contourf(X,Y,Z,100,cmap='Greens')
plt.colorbar()
ax.set_title('Bz')
ax.set_xlabel('Freq (GHz)')
ax.set_ylabel('B (uT)')

ax = fig.add_subplot(2,3,3)
X = freq/1e9
##Y = spectraxy[:,0]*1e6
Y = spectraxz[:,0]*1e6
X,Y = np.meshgrid(X,Y)
##Z = spectraxy[:,1:]
Z = spectraxz[:,1:]
plt.contourf(X,Y,Z,100,cmap='Blues')
plt.colorbar()
ax.set_title('Bxz')
ax.set_xlabel('Freq (GHz)')
ax.set_ylabel('B (uT)')

ax = fig.add_subplot(2,3,4)
X = freq/1e9
##Y = spectraxy[:,0]*1e6
Y = spectraxyz[:,0]*1e6
X,Y = np.meshgrid(X,Y)
##Z = spectraxy[:,1:]
Z = spectraxyz[:,1:]
plt.contourf(X,Y,Z,100,cmap='Greys')
plt.colorbar()
ax.set_title('Bxyz')
ax.set_xlabel('Freq (GHz)')
ax.set_ylabel('B (uT)')

ax = fig.add_subplot(2,3,5)
ax.plot(ewx[:,1]/1e6,spectrax[:,0]*1e6,'r-',\
        ewz[:,1]/1e6,spectraz[:,0]*1e6,'g-',\
        ewxz[:,1]/1e6,spectraxz[:,0]*1e6,'b-',\
        ewxyz[:,1]/1e6,spectraxyz[:,0]*1e6,'k-')
ax.set_title('Splitting')
ax.set_xlabel('Freq (MHz)')
ax.set_ylabel('B (uT)')

ax = fig.add_subplot(2,3,6)
ax.plot(ewx[:,2]/1e6,spectrax[:,0]*1e6,'r-',\
        ewz[:,2]/1e6,spectraz[:,0]*1e6,'g-',\
        ewxz[:,2]/1e6,spectraxz[:,0]*1e6,'b-',\
        ewxyz[:,2]/1e6,spectraxyz[:,0]*1e6,'k-')
ax.set_title('Width')
ax.set_xlabel('Freq (MHz)')
ax.set_ylabel('B (uT)')

plt.subplots_adjust(wspace=.2,hspace=1)
plt.draw()
plt.show()

