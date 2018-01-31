import numpy as np
import matplotlib.pyplot as plt

def Bdom(doms,z,n):
    xLength = 0.001
    yLength = 0.001
    zLength = 0.001

    xstep = xLength/doms
    xArray = np.arange(-xLength/2+xstep/2,xLength/2,xstep)
##    print 'xstep = ' + str(xstep)
##    print 'xArray = ' + str(xArray)

    dV = xstep*yLength*zLength
    mu = 1.71 * (9.274*(10**(-24))) # J / (T atom)
    rho = 8900000. # kg / m3 - real value: 8900000 kg/m3
    molmass = 58.9332 # g / mol
    Navo = 6.022*(10**23) # atoms / mol
    M = (mu * Navo * rho) / molmass
    m = np.array([0,0,M*dV])
##    print 'mz = ' + str(mz)

    def Bx(doms,n):
        mu0 = 4*np.pi*(10**(-7)) # T m / A
        r = np.array([-xArray[n],0,z])
        rmag = np.linalg.norm(r)
        rhat = r/rmag
        mn = (m*((-1)**(n+1)))
##        print 'r = ' + str(r)
##        print 'm = ' + str(mn)
##        print 'm dot rhat = ' + str(np.inner(mn,rhat))
        return (mu0/(4.*np.pi))\
               *(3.*np.inner(mn,rhat)*rhat-mn)\
               /(rmag**3)

    return Bx(doms,n)

doms = np.zeros((200))
for i in range(0,len(doms),1):
    doms[i] = 10*(i+1)+540
print doms
##for dom in [1,2,4,8,16,32,64,128,256,512,1024,2048,4096]:
##doms = np.array([2,4])
##Bsum = np.zeros((len(doms),4))
##for i in range(0,len(doms),1):
##    for j in range(0,int(doms[i]),1):
##        Bsum[i,0] = doms[i]
##        Bsum[i,1:] = Bsum[i,1:] + Bdom(doms[i],-0.001,j)
Bsum = np.zeros((len(doms),5))
for i in range(0,len(doms),1):
    for j in range(0,int(doms[i]),1):
        Bsum[i,0] = doms[i]
        Bsum[i,1] = .001/doms[i]
        Bsum[i,2:] = Bsum[i,2:] + Bdom(doms[i],-0.001,j)




print Bsum


np.savetxt('domain-Bx-By-Bz-array.txt',\
           Bsum, delimiter = ', ')











