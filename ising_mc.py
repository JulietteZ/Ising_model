#%%
import numpy as np
import matplotlib.pyplot as plt
import math

#%%
# Build the system
def hot_start():
    lattice = np.random.randint(0,2,(ns,ns))
    lattice[lattice==0] = -1
    return lattice

def cold_start():
    lattice = np.ones((ns,ns))
    return lattice

# Periodic boundary conditions
def bc(i):
    if i > ns-1:
        return 0
    if i < 0:
        return ns-1
    else:
        return i

# Measure magnetization
def mag(lattice):
    m = 0.0
    for j in range(0,ns):
        for k in range(0,ns):
            m += lattice[j,k]
    return m/(ns*ns)

def afm(lattice, coordinate_n=4):
    m = 0.
    for j in range(0,ns):
        for k in range(0,ns):
            m += energy(lattice, j, k)/coordinate_n
    return m/(ns*ns)   

# Calculate internal energy
def energy(lattice, N, M):
    return -1 * lattice[N,M] * (lattice[bc(N-1), M] \
                       + lattice[bc(N+1), M] \
                       + lattice[N, bc(M-1)] \
                       + lattice[N, bc(M+1)])

def sum_nn(lattice, N, M):
    return (lattice[bc(N-1), M] + lattice[bc(N+1), M] + lattice[N, \
                    bc(M-1)] + lattice[N, bc(M+1)])


# The Main monte carlo loop
def update(beta):
    lattice = hot_start()
    for step in enumerate(range(ns*ns)):
        j = np.random.randint(0,ns)
        k = np.random.randint(0,ns)
        
        E = -2. * energy(lattice, j, k)
        
        if E <= 0.:
            lattice[j,k] *= -1
        elif np.exp(-beta*E) > np.random.rand():
            lattice[j,k] *= -1
            
def sweep(lattice, beta):
    acc = 0
    for j in range(0,ns):
        for k in range(0,ns):
            new_spin = -lattice[j,k]
            dE = 1*(new_spin-lattice[j,k])*sum_nn(lattice, j, k)
            if dE <= 0.:
                lattice[j,k] = new_spin
                acc += 1
            elif np.exp(-beta*dE) > np.random.rand():
                lattice[j,k] = new_spin
                acc += 1
    accept = (1.*acc)/(ns*ns)

#%%
# Main
if __name__ == "__main__":
    ns = 20
    ninit = 300
    nsweeps = 1000

    for size in range(1):
        mav = np.zeros(40)
        erg = np.zeros(nsweeps)
        heat = np.zeros(40)
        beta = np.ones(40)
        for b in range(40):
            beta[b] = b/100+0.2
            lattice = hot_start()
            for n in range(ninit):
                sweep(lattice,beta[b])
            mav[b] = 0
            erg[0] = 0
            mlist = np.ones(nsweeps)
            for n in range(nsweeps):
                sweep(lattice,beta[b])
                m = afm(lattice)
                mav[b] += m
                erg[n] = afm(lattice) * 2.
                mlist[n] = m
            mav[b] = mav[b]/nsweeps
            heat[b] = np.sum(pow(erg-np.mean(erg),2))/nsweeps
        
        plt.plot(mlist)
        plt.title("afm square lattice size: %i * %i, beta: %f" %(ns,ns,beta[-1]))
        plt.xlabel("nsweeps")
        plt.ylabel("average m")
        plt.show()

        plt.plot(beta,heat)
        plt.title("afm square lattice size: %i * %i" %(ns,ns))
        plt.xlabel(r'$\beta$')
        plt.ylabel("specific heat")
        plt.show()

        plt.plot(beta,np.abs(mav))
        plt.title("afm square lattice size: %i * %i" %(ns,ns))
        plt.xlabel(r'$\beta$')
        plt.ylabel("average m")
        plt.show()
        ns = ns*2

# %%
