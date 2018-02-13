import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import datetime

from NVeigenvalues import eigenvalues
from lor8 import lor8
from DEwidthZeroFieldFit import de

#Load arrays [# of doms, dom size, Bx,By,Bz]
dBxyzArray = np.loadtxt('domain-Bx-By-Bz-array.txt',\
             delimiter=', ',unpack=False)
sim1600um = np.loadtxt('z_-1600um.txt',\
             delimiter=', ',unpack=False)
sim1400um = np.loadtxt('z_-1400um.txt',\
             delimiter=', ',unpack=False)
sim1200um = np.loadtxt('z_-1200um.txt',\
             delimiter=', ',unpack=False)
sim800um = np.loadtxt('z_-800um.txt',\
             delimiter=', ',unpack=False)
sim600um = np.loadtxt('z_-600um.txt',\
             delimiter=', ',unpack=False)
sim400um = np.loadtxt('z_-400um.txt',\
             delimiter=', ',unpack=False)
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
# Fit the Co data for adding splitting and width
# to plots below
zfArraytmp = np.array([2.87e9,3e6,2e6,0,0,0,0,0])
CoFitParams = np.copy(de(freqCo,signalCo,zfArraytmp))
NVFitParams = np.copy(de(freqNV,signalNV,zfArraytmp))







# Take D, E, and width from Co zf data for sim
zfArray = np.zeros(8)
zfArray[0:3] = np.copy(NVFitParams[0:3])
# Add EV columns to Bxyz arrays
dBxyzArray = np.column_stack((dBxyzArray,
             np.zeros((len(dBxyzArray),7))))
sim1600um = np.column_stack((sim1600um,
             np.zeros((len(sim1600um),7))))
sim1400um = np.column_stack((sim1400um,
             np.zeros((len(sim1400um),7))))
sim1200um = np.column_stack((sim1200um,
             np.zeros((len(sim1200um),7))))
sim800um = np.column_stack((sim800um,
             np.zeros((len(sim800um),7))))
sim600um = np.column_stack((sim600um,
             np.zeros((len(sim600um),7))))
sim400um = np.column_stack((sim400um,
             np.zeros((len(sim400um),7))))
# Set Bz = 0
dBxyzArray[:,4] = 0
sim1600um[:,4] = 0
sim1400um[:,4] = 0
sim1200um[:,4] = 0
sim800um[:,4] = 0
sim600um[:,4] = 0
sim400um[:,4] = 0
##print dBxyzArray


###################################################################
###################################################################
# This is an aside...
# The interesting thing is that when dealing with Bz, you get
# two peaks at large fields... however, when dealing with Bx,
# you get more peaks at high fields... this must have something
# to do with the (Sx^2+Sy^2) term in the Hamiltonian???
##Bz = np.array([[0,0,1e-2],[1e-2,0,0]])
##zf = np.array([2.87e9,3e6,2e6,0])
##EV = np.zeros((len(Bz),8))
##for i in range(0,len(Bz),1):
##    EV[i,:] = eigenvalues(zf,np.array(Bz[i]))
####    EV[i,:] = eigenvalues(zf,np.array([0,0,Bz[i]]))
####print EV
##freq = np.arange(2.5e9,3.25e9,1e6)
##amp = np.array([1e7,1e7,1e7,1e7,1e7,1e7,1e7,1e7])
##plt.plot(freq,lor8(freq,zf,amp,EV[0]),\
##         freq,lor8(freq,zf,amp,EV[1]))
##plt.title('blue = large Bz, orange = large Bx',fontsize=15)
###################################################################
###################################################################
###################################################################


#calculate eigenvalues for each number of domains
for i in range(0,len(dBxyzArray),1):
    simulatedEV = eigenvalues(zfArray[0:3],dBxyzArray[i,1:4])
    dBxyzArray[i,4:12] = simulatedEV
    ev1600um = eigenvalues(zfArray[0:3],sim1600um[i,1:4])
    sim1600um[i,4:12] = ev1600um
    ev1400um = eigenvalues(zfArray[0:3],sim1400um[i,1:4])
    sim1400um[i,4:12] = ev1400um
    ev1200um = eigenvalues(zfArray[0:3],sim1200um[i,1:4])
    sim1200um[i,4:12] = ev1200um
    ev800um = eigenvalues(zfArray[0:3],sim800um[i,1:4])
    sim800um[i,4:12] = ev800um
    ev600um = eigenvalues(zfArray[0:3],sim600um[i,1:4])
    sim600um[i,4:12] = ev600um
    ev400um = eigenvalues(zfArray[0:3],sim400um[i,1:4])
    sim400um[i,4:12] = ev400um
# Isolate eigenvalue arrays
evArray = np.delete(dBxyzArray,(0,1,2,3),axis=1)
print dBxyzArray[0]
print evArray[0]
print sim800um[0]
ev1600um = np.delete(sim1600um,(0,1,2,3),axis=1)
ev1400um = np.delete(sim1400um,(0,1,2,3),axis=1)
ev1200um = np.delete(sim1200um,(0,1,2,3),axis=1)
ev800um = np.delete(sim800um,(0,1,2,3),axis=1)
ev600um = np.delete(sim600um,(0,1,2,3),axis=1)
ev400um = np.delete(sim400um,(0,1,2,3),axis=1)
print ev800um[0]

###################################################################
###################################################################
###################################################################
# Side note: spectra from 20G helmholtz coil
##helmholtzB = np.array([0.0018,0.0018,0.0018]) # 2mT = 20G
##helmholtzEV = eigenvalues(zfArray[0:3],helmholtzB)
##helmholtzFreq = np.arange(2.57e9,3.17e9,1e5)
##helmholtzAmp = np.array([1e7,1e7,1e7,1e7,1e7,1e7,1e7,1e7])
##helmholtzSpectra = lor8(helmholtzFreq,zfArray,\
##                        helmholtzAmp,helmholtzEV)
##plt.plot(helmholtzFreq,helmholtzSpectra)
##plt.show()
###################################################################
###################################################################
###################################################################


# Generate spectra for each set of domains
freq = np.arange(2.845e9,2.895e9,1e5)
ampArray = np.array([1e7,1e7,1e7,1e7,1e7,1e7,1e7,1e7])
spectraArray = np.zeros((len(dBxyzArray),len(freq)+1))
spectra1600um = np.zeros((len(sim1600um),len(freq)+1))
spectra1400um = np.zeros((len(sim1400um),len(freq)+1))
spectra1200um = np.zeros((len(sim1200um),len(freq)+1))
spectra800um = np.zeros((len(sim800um),len(freq)+1))
spectra600um = np.zeros((len(sim600um),len(freq)+1))
spectra400um = np.zeros((len(sim400um),len(freq)+1))
print 
for i in range(0,len(dBxyzArray),1):
    spectraArray[i,0] = dBxyzArray[i,1]*1e6
    spectraArray[i,1:] = lor8(freq,zfArray,\
                             ampArray,evArray[i,:])
    spectraArray[i,1:] = spectraArray[i,1:]/\
                         np.amax(spectraArray[i,1:])
    
    spectra1600um[i,0] = sim1600um[i,1]*1e6
    spectra1600um[i,1:] = lor8(freq,zfArray,\
                             ampArray,ev1600um[i,:])
    spectra1600um[i,1:] = spectra1600um[i,1:]/\
                         np.amax(spectra1600um[i,1:])

    spectra1400um[i,0] = sim1400um[i,1]*1e6
    spectra1400um[i,1:] = lor8(freq,zfArray,\
                             ampArray,ev1400um[i,:])
    spectra1400um[i,1:] = spectra1400um[i,1:]/\
                         np.amax(spectra1400um[i,1:])

    spectra1200um[i,0] = sim1200um[i,1]*1e6
    spectra1200um[i,1:] = lor8(freq,zfArray,\
                             ampArray,ev1200um[i,:])
    spectra1200um[i,1:] = spectra1200um[i,1:]/\
                         np.amax(spectra1200um[i,1:])

    spectra800um[i,0] = sim800um[i,1]*1e6
    spectra800um[i,1:] = lor8(freq,zfArray,\
                             ampArray,ev800um[i,:])
    spectra800um[i,1:] = spectra800um[i,1:]/\
                         np.amax(spectra800um[i,1:])

    spectra600um[i,0] = sim600um[i,1]*1e6
    spectra600um[i,1:] = lor8(freq,zfArray,\
                             ampArray,ev600um[i,:])
    spectra600um[i,1:] = spectra600um[i,1:]/\
                         np.amax(spectra600um[i,1:])

    spectra400um[i,0] = sim400um[i,1]*1e6
    spectra400um[i,1:] = lor8(freq,zfArray,\
                             ampArray,ev400um[i,:])
    spectra400um[i,1:] = spectra800um[i,1:]/\
                         np.amax(spectra400um[i,1:])

##print spectraArray
##print dBxyzArray
EArray = np.zeros(len(spectraArray))
EArrayErr = np.zeros(len(spectraArray))
widthArray = np.zeros(len(spectraArray))
widthArrayErr = np.zeros(len(spectraArray))
E1600um = np.zeros(len(spectra1600um))
width1600um = np.zeros(len(spectra1600um))
E1400um = np.zeros(len(spectra1400um))
width1400um = np.zeros(len(spectra1400um))
E1200um = np.zeros(len(spectra1200um))
width1200um = np.zeros(len(spectra1200um))
E800um = np.zeros(len(spectra800um))
width800um = np.zeros(len(spectra800um))
E600um = np.zeros(len(spectra800um))
width600um = np.zeros(len(spectra800um))
E400um = np.zeros(len(spectra800um))
width400um = np.zeros(len(spectra800um))
for i in range(0,len(spectraArray),1):
    zfTmp = de(freq,spectraArray[i,1:],\
                zfArray)
    EArray[i] = zfTmp[1]
    EArrayErr[i] = zfTmp[5] # error on E
    widthArray[i] = zfTmp[2]
    widthArrayErr[i] = zfTmp[6] # Error on width
#---------------------------------------------
    zf1600um = de(freq,spectra1600um[i,1:],\
                zfArray)
    E1600um[i] = zf1600um[1]
    width1600um[i] = zf1600um[2]
#---------------------------------------------
    zf1400um = de(freq,spectra1400um[i,1:],\
                zfArray)
    E1400um[i] = zf1400um[1]
    width1400um[i] = zf1400um[2]
#---------------------------------------------
    zf1200um = de(freq,spectra1200um[i,1:],\
                zfArray)
    E1200um[i] = zf1200um[1]
    width1200um[i] = zf1200um[2]
#---------------------------------------------
##    zf800um = de(freq,spectra800um[i,1:],\
##                zfArray)
##    E800um[i] = zf800um[1]
##    width800um[i] = zf800um[2]
##print E800um[0]
print np.column_stack([dBxyzArray[:,1]*1e6,\
                       EArray,EArrayErr,\
                       widthArray,widthArrayErr])





def domainSizeVsBx(dBxyzArray,freq,spectraArray):
    domain_size = dBxyzArray[:,1]*1e6 # units um
    Bx = dBxyzArray[:,2]*1e6 # units: uT
    Bz = dBxyzArray[:,4]*1e6 # units: uT
    Bx1600um = sim1600um[:,2]*1e6
    Bx1400um = sim1400um[:,2]*1e6
    Bx1200um = sim1200um[:,2]*1e6
    fig = plt.figure(figsize=plt.figaspect(1.))

    ax = fig.add_subplot(221)
    ax.plot(domain_size,Bx,'r-',\
            domain_size,Bx1600um,'r:',\
            domain_size,Bx1400um,'r-.',\
            domain_size,Bx1200um,'r--')
    ax.set_title('Bx as a function of domain size (z = 1mm, 1.2mm, 1.4mm, 1.6mm)')
    ax.set_xlabel('Domain size (um)')
    ax.set_ylabel('Bx (uT)')
    ax.grid(True)
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
#--------------------------------------------------------
    ax = fig.add_subplot(223)
    horiz_line_data = np.array([1.75 for i in xrange(len(EArray))])
    ax.plot(EArray,dBxyzArray[:,1]*1e6,'r-',\
            widthArray,dBxyzArray[:,1]*1e6,'b-',\
            E1200um,sim1200um[:,1]*1e6,'r--',\
            width1200um,sim1200um[:,1]*1e6,'b--',\
            E1400um,sim1200um[:,1]*1e6,'r-.',\
            width1400um,sim1400um[:,1]*1e6,'b-.',\
            E1600um,sim1200um[:,1]*1e6,'r:',\
            width1600um,sim1600um[:,1]*1e6,'b:',\
            EArray,horiz_line_data,'g--')
    ax.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    ax.set_title('Splitting (red) and width (blue) as a function of domain size')
    ax.set_xlabel('Freq (Hz)')
    ax.set_ylabel('Domain size (um)')
    ax.grid(True)
#--------------------------------------------------------
    ax = fig.add_subplot(224)
    X = freq
    Y = dBxyzArray[:,1]*1e6
##    Y = sim600um[:,1]*1e6
    X, Y = np.meshgrid(X, Y)
    Z = spectraArray[:,1:]
##    Z = spectra600um[:,1:]
    plt.contourf(X,Y,Z,100,cmap='RdGy')
    plt.colorbar()    
    ax.set_title('NV spectra as a function of domain size (z=-1mm)')
##    ax.set_title('NV spectra as a function of domain size (z=-0.6mm)')
    ax.set_ylabel('Domain size (um)')
    ax.set_xlabel('Frequency (Hz)')
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











