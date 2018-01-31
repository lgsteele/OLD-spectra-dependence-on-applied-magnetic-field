import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import datetime

from NVeigenvalues import eigenvalues
from lor8 import lor8
from DEwidthZeroFieldFit import de

#Load array [dom,Bx,By,Bz]
dBxyzArray = np.loadtxt('domain-Bx-By-Bz-array.txt',\
             delimiter=', ',unpack=False)
##print dBxyzArray

# Upload Cobalt data
freqCo,signalCo = np.loadtxt('Co-sum.txt',delimiter=', ',\
                             skiprows=1,unpack=True)
signalCo = signalCo+450
freqCofit,signalCofit = np.loadtxt('Co-fitArray.txt',delimiter=', ',\
                             skiprows=1,unpack=True)
signalCofit = signalCofit+450
freqNV,signalNV = np.loadtxt('NV-after-sum.txt',delimiter=', ',\
                             skiprows=1,unpack=True)
freqNVfit,signalNVfit = np.loadtxt('NV-fitArray.txt',delimiter=', ',\
                             skiprows=1,unpack=True)
# Fit the spectra
zfArraytmp = np.array([2.87e9,3e6,2e6,0,0,0,0,0])
CoFitParams = np.copy(de(freqCo,signalCo,zfArraytmp))
NVFitParams = np.copy(de(freqNV,signalNV,zfArraytmp))







# Plot peak splitting and width (D & E) in low field limit
#add eigenvalues to array
zfArray = np.zeros(8)
zfArray[0:3] = np.copy(NVFitParams[0:3])
print zfArray
##zfArray = np.array([2.87e9,3e6,2e6,0,0,0,0,0])
dBxyzArray = np.column_stack((dBxyzArray,
             np.zeros((len(dBxyzArray),8))))
dBxyzArray[:,4] = 0
##print dBxyzArray


###################################################################
###################################################################
# This is an aside...
# The interesting thing is that when dealing with Bz, you get
# two peaks at large fields... however, when dealing with Bx,
# you get more peaks at high fields... this must have something
# to do with the (Sx^2+Sy^2) term in the Hamiltonian???
Bz = np.array([[0,0,1e-2],[1e-2,0,0]])
zf = np.array([2.87e9,3e6,2e6,0])
EV = np.zeros((len(Bz),8))
for i in range(0,len(Bz),1):
    EV[i,:] = eigenvalues(zf,np.array(Bz[i]))
##    EV[i,:] = eigenvalues(zf,np.array([0,0,Bz[i]]))
print EV
freq = np.arange(1.87e9,3.87e9,1e6)
amp = np.array([1e7,1e7,1e7,1e7,1e7,1e7,1e7,1e7])
plt.plot(freq,lor8(freq,zf,amp,EV[0]),\
         freq,lor8(freq,zf,amp,EV[1]))
###################################################################
###################################################################
###################################################################


#calculate eigenvalues for each number of domains
for i in range(0,len(dBxyzArray),1):
    simulatedEV = eigenvalues(zfArray[0:3],dBxyzArray[i,1:4])
    dBxyzArray[i,4:12] = simulatedEV
# Isolate eigenvalue array
evArray = np.delete(dBxyzArray,(0,1,2,3),axis=1)
# Generate spectra for each set of domains
freq = np.arange(2.845e9,2.895e9,1e5)
ampArray = np.array([1e7,1e7,1e7,1e7,1e7,1e7,1e7,1e7])
spectraArray = np.zeros((len(dBxyzArray),len(freq)+1))
for i in range(0,len(dBxyzArray),1):
    spectraArray[i,0] = dBxyzArray[i,1]*1e6
    spectraArray[i,1:] = lor8(freq,zfArray,\
                             ampArray,evArray[i,:])
    spectraArray[i,1:] = spectraArray[i,1:]/\
                         np.amax(spectraArray[i,1:])
##print spectraArray
##print dBxyzArray
EArray = np.zeros(len(spectraArray))
widthArray = np.zeros(len(spectraArray))
for i in range(0,len(spectraArray),1):
    zfTmp = de(freq,spectraArray[i,1:],\
                zfArray)
    EArray[i] = zfTmp[1]
    widthArray[i] = zfTmp[2]







def domainSizeVsBx(dBxyzArray,freq,spectraArray):
    domain_size = dBxyzArray[:,1]*1e6 # units um
    Bx = dBxyzArray[:,2]*1e6 # units: uT
    Bz = dBxyzArray[:,4]*1e6 # units: uT

    fig = plt.figure(figsize=plt.figaspect(1.))

    ax = fig.add_subplot(221)
    ax.plot(domain_size,Bx)
    ax.set_title('Bx as a function of domain size')
    ax.set_xlabel('Domain size (um)')
    ax.set_ylabel('Bx (uT)')
    ax.grid(True)

##    ax = fig.add_subplot(222)
##    ax.plot(domain_size,Bx)
##    ax.set_title('Bx as a function of domain size (log)')
##    ax.set_xscale('log')
##    ax.set_xlabel('Domain size (um)')
##    ax.set_ylabel('Bx (uT)')
##    ax.grid(True)
    ax = fig.add_subplot(223)
    horiz_line_data = np.array([1.75 for i in xrange(len(EArray))])
    ax.plot(EArray,dBxyzArray[:,1]*1e6,'r-',\
            widthArray,dBxyzArray[:,1]*1e6,'b-',\
            EArray,horiz_line_data,'g--')
    ax.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    ax.set_title('Splitting (red) and width (blue) as a function of domain size')
    ax.set_xlabel('Freq (Hz)')
    ax.set_ylabel('Domain size (um)')
    ax.grid(True)
##    ax.axis([2,6,.5,3])


    ax = fig.add_subplot(224)
    X = freq
##    Y = spectraArray[:,0]
    Y = dBxyzArray[:,1]*1e6
    X, Y = np.meshgrid(X, Y)
    Z = spectraArray[:,1:]
    plt.contourf(X,Y,Z,100,cmap='RdGy')
    plt.colorbar()    
    ax.set_title('NV spectra as a function of domain size')
    ax.set_ylabel('Domain size (um)')
    ax.set_xlabel('Frequency (Hz)')
#------------------------
    ax = fig.add_subplot(222)
    ax.set_title('NV spectra with and without 1m3 Co sample')
    ax.set_xlabel('frequency (Hz)')
    ax.set_ylabel('signal (AU)')
    ax.text(2.88e9,200,'w/o Co \nsplitting=%sMHz \nwidth=%sMHz'%\
            (round(NVFitParams[1]/1e6,2),\
             round(NVFitParams[2]/1e6,2)))
    ax.text(2.88e9,700,'w/ Co \nsplitting=%sMHz \nwidth=%sMHz'%\
            (round(CoFitParams[1]/1e6,2),\
             round(CoFitParams[2]/1e6,2)))
##    ax.text(2.804e9,800,'Hz = (E_Co - E_NV)/(2*gamma)')
##    ax.text(2.804e9,700,'Hz = 25.725 +/- 10.186 uT')
##    ax.text(2.804e9,600,'dHz\' = dHz + (w-w0)/2')
    ax.plot(freqCo,signalCo,'b--',\
             freqCofit,signalCofit,'r-',\
             freqNV,signalNV,'b--',\
             freqNVfit,signalNVfit,'r-')
    plt.yticks([])
##    ax.plot(domain_size,Bx)
##    ax.set_yscale('log')
##    ax.set_xlabel('Domain size (um)')
##    ax.set_ylabel('Bx (uT)')
##    ax.grid(True)
##    plt.tight_layout()
    plt.subplots_adjust(wspace=.2,hspace=.3)
    plt.draw()
    plt.show()
##    plt.plot(domain_size,Bx,'ro')
##    plt.xscale('log')
##    plt.yscale('log')
##    plt.show()
##    print domain_size
domainSizeVsBx(dBxyzArray,freq,spectraArray)











