import numpy as np
import matplotlib.pyplot as plt
import os, glob, sys

filename = glob.glob('*.txt')
m = int(len(filename))
print m
def data():
    """
    Import data files and put the eigenvectors into arrays.
    """
    rho = np.loadtxt(filename[-1])[:,0]
    N = len(rho)
    Eigvecs = np.zeros((N, m))
    twoElectrons = np.zeros((N, m/2))
    for ind, name in enumerate(filename):
        data = np.loadtxt(name)
        Eigvecs[:,ind] = data[:,1]

    return rho, Eigvecs

def plotEigenvectors():
    """
    Plot the ground stages
    """
    rho, Eigvecs = data()
    print np.shape(Eigvecs)
    color = ['b', 'r', 'g', 'm', 'b', 'r', 'g', 'm']
    omega_r = [0.01, 0.5, 1.0, 5.0, 0.01, 0.5, 1.0, 5.0]

    for i in range(m):
        if i >= m/2:
            plt.figure('2 electrons with repulsion')
            plt.plot(rho, abs(Eigvecs[:,i])**2, c=color[i], label=r'$\omega_r$=%d'%(omega_r[i]))
            plt.xlabel(r'$\rho$', size=14)
            plt.ylabel(r'$|\psi(\rho)|^2$', size=14)
            plt.legend(loc=1, fontsize=12)
            plt.grid('on')
            #plt.savefig('2el_repulsion.png')

        else:
            plt.figure('2 electrons, No repulsion')
            plt.plot(rho, abs(Eigvecs[:,i])**2,c=color[i], label=r'$\omega_r$=%d'%(omega_r[i]))
            plt.xlabel(r'$\rho$', size=14)
            plt.ylabel(r'$|\psi(\rho)|^2$', size=14)
            plt.legend(loc=1, fontsize=12)
            plt.grid('on')
            #plt.savefig('2el_norepulsion.png')


d = data()
plot = plotEigenvectors()
plt.show()
sys.exit()
try:
    filename = sys.argv[1]#'sch1_n20.txt'#glob.glob('sch1_n20.txt')
except:
    print 'No input file. Type filename to be loaded'
    sys.exit()
print


#sys.exit()
def Plotdata(filename):
    data = np.loadtxt(filename)
    N = len(data)
    print np.shape(data)
    rho = data[:,0]
    Eigenvectors = np.array(data[:,1:])
    print len(Eigenvectors[0])

    # Plotting:
    plt.figure('1d Scroedinger')
    color = ['b', 'r', 'g']
    for i in range(len(Eigenvectors[0])):
        if i == 0:
            plt.plot(rho, abs(Eigenvectors[:,i])**2, c=color[0], label='Ground state')
        else:
            plt.plot(rho, abs(Eigenvectors[:,i])**2, c=color[i], label='n = %d'%(i))

    plt.xlabel(r'$\rho$', size=14)
    plt.ylabel(r'$|\psi(\rho)|^2$', size=14)
    plt.legend(loc=1, fontsize=12)
    plt.grid('on')
    plt.savefig('Schr_1d_n%d.png'%(len(rho)))



#o= Plotdata(filename)
plt.show()
