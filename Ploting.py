import numpy as np
import matplotlib.pyplot as plt
import os, glob, sys



filenames = glob.glob('*.txt')

repulsion = []; norep = []; states = []
for ind, name in enumerate(filenames):
    if name == 'sch2omega_r%d.txt'%(ind+1-8):
        #print name
        repulsion.append(name)
    if name == 'norepulsionomega_r%d.txt'%(ind+1):
        #print name
        norep.append(name)
    if name == 'sch1_%d.txt'%(ind+1-4):
        #print name
        states.append(name)

#sys.exit()
def data2electrons(filename):
    """
    Import data files and put the eigenvectors into arrays.
    """
    print filename
    rho = np.loadtxt(filename[-1])[:,0]
    N = len(rho); m = len(filename)
    print N, m
    Eigvecs = np.zeros((N, m))
    #repulsion = []; norep = []
    for ind, name in enumerate(filename):
        data = np.loadtxt(name)
        #print np.shape(data[:,1]), np.shape(Eigvecs[:,ind])
        Eigvecs[:,ind] = data[:,1]

    return rho, Eigvecs

def dataStates(filename):
    """
    Load .txt files for one electron into arrays.
    """
    rho = np.loadtxt(filename[-1])[:,0]
    N = len(rho); m = len(filename)
    Eigvec1 = np.zeros((N, m-1))
    Eigvec2 = np.zeros((N, m-1))
    Eigvec3 = np.zeros((N, m-1))
    Eigvec4 = np.zeros((N, m-1))
    Eigvecs = np.zeros((m, N, m-1))#np.array([Eigvec1,Eigvec2, Eigvec3, Eigvec4])
    print np.shape(Eigvecs)
    for i in range(m):
        Eigvecs[i,:,:] = np.loadtxt(filename[i])[:,1:]

    return rho, Eigvecs


def plotStates(filename, name):
    rho, Eigvecs = dataStates(filename)
    omega_r = [0.01, 0.5, 1.0,5.0]
    col = ['b', 'r', 'g']; wr = ['001', '05', '1','5']
    for w in range(len(omega_r)):
        plt.figure('States for %g'%(omega_r[w]))
        for i in range(len(Eigvecs[1,1,:])):
            plt.plot(rho, abs(Eigvecs[w,:,i])**2, c=col[i], label='n=%d'%i)

        plt.legend(loc=1, fontsize=12)
        plt.xlabel(r'$\rho$', size=14)
        plt.ylabel(r'$|\psi(\rho)|^2$', size=14)
        plt.tight_layout()
        plt.savefig('Plot_pro2/States_%s.png'%(wr[w]))


def plotEigenvectors(filename, name):
    """
    Plot the ground stages
    """
    rho, Eigvecs = data2electrons(filename)
    print np.shape(Eigvecs)
    color = ['b', 'r', 'g', 'm', 'b', 'r', 'g', 'm']
    omega_r = [0.01, 0.5, 1.0, 5.0, 0.01, 0.5, 1.0, 5.0]
    #deigvec = np.gradient(abs(Eigvecs[:,4])**2)
    m = len(filename)
    for i in range(m):
        #if i >= m/2:
        plt.figure('%s s_normalised'%name)
        plt.plot(rho, (np.absolute(Eigvecs[:,i])**2)/(np.max(np.absolute(Eigvecs[:,i])**2)), \
                c=color[i], label=r'$\omega_r$=%g'%(omega_r[i]))
        plt.xlabel(r'$\rho$', size=14)
        plt.ylabel(r'$|\psi(\rho)|^2$', size=14)
        plt.legend(loc=9, fontsize=12)
        plt.grid('on')
        plt.tight_layout()
        plt.savefig('Plot_pro2/2el_%s_normalised.png'%(name))
        plt.savefig('Plot_pro2/2el_%s_normalised_rho10.png'%(name))

        plt.figure('%s'%name)
        plt.plot(rho, np.absolute(Eigvecs[:,i])**2,c=color[i], label=r'$\omega_r$=%g'%(omega_r[i]))
        plt.xlabel(r'$\rho$', size=14)
        plt.ylabel(r'$|\psi(\rho)|^2$', size=14)
        plt.legend(loc=9, fontsize=12)
        plt.grid('on')
        plt.tight_layout()
        #plt.savefig('Plot_pro2/2el_%s.png'%name)
        plt.savefig('Plot_pro2/2el_%s_rho10.png'%name)


#d = data()
plot1 = plotEigenvectors(repulsion, 'repulsion')
#p2 = plotEigenvectors(norep, 'no_repulsion')
#p3 = plotStates(states, 'states')

plt.show()
