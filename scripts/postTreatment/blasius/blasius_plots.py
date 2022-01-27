#!/usr/bin/python3

#-------------------------------------------------------------------------------------------------------------------------------
#                                                           OBJECTIVE                                                           
#-------------------------------------------------------------------------------------------------------------------------------
# This script (2/2) is for the post-treatment of a subsonic flow over a flat flate, i.e. Blasius' test case.
# It aims to compare steady dimensionless normal/tangential velocities with the analytical solution.
# The analytical solution is based on the book:
#  'I do like CFD', vol. 1, Katate Masatsuka
# The comparison is done through a convergence study with L2 norms.

#-------------------------------------------------------------------------------------------------------------------------------
#                                                            IMPORTS                                                            
#-------------------------------------------------------------------------------------------------------------------------------
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import math

#-------------------------------------------------------------------------------------------------------------------------------
#                                           PARAMETERS TO DEFINE / INITIAL CONDITIONS                                           
#-------------------------------------------------------------------------------------------------------------------------------
# Observables
# -----------
files = ['res_mesh_200x90.txt', 'res_mesh_400x105.txt', 'res_mesh_800x120.txt', 'res_mesh_1600x135.txt']
    # Filenames of post-treated results 
    # (simulation w/ corresponding analytical sol.)
file_exact = 'exact.txt' # Analytical sol. for plot (continuous)

#-------------------------------------------------------------------------------------------------------------------------------
#                                          READING RESULTS FROM Blasius_computation.py                                                
#-------------------------------------------------------------------------------------------------------------------------------
# Post-treated results
# --------------------
Nmesh = len(files) # Number of test cases
y = []             # Array of vertical positions
ustar = []         # Array of dimensionless normal velocity
vstar = []         # Array of dimensionless tangential velocity
eta = []           # Array of Blasius dimensionless variable
u_bls = []         # Array of dimensionless analytical normal velocity
v_bls = []         # Array of dimensionless analytical tangential velocity

L2_analytical = np.zeros(Nmesh) # L2 norm array build upon analytical sol.
dy = np.zeros(Nmesh)            # Array of vertical size of 1st cell next to the wall
dy_str = []                     # Array of vertical size for plot

# Read data of each mesh
# Column names
# 0. y
# 1. eta
# 2. ustar
# 3. vstar
# 4. u_bls
# 5. v_bls
# 6. L2_analytical

count = 0
for f in files:
    y.append(np.loadtxt(f, delimiter=' ', skiprows=1)[:,0])
    eta.append(np.loadtxt(f, delimiter=' ', skiprows=1)[:,1])
    ustar.append(np.loadtxt(f, delimiter=' ', skiprows=1)[:,2])
    vstar.append(np.loadtxt(f, delimiter=' ', skiprows=1)[:,3])
    u_bls.append(np.loadtxt(f, delimiter=' ', skiprows=1)[:,4])
    v_bls.append(np.loadtxt(f, delimiter=' ', skiprows=1)[:,5])
    L2_analytical[count] = (np.loadtxt(f, delimiter=' ', skiprows=1)[0,6])
    count += 1
    
# Creates array of dy for x-axis of L2 plot
dy = np.zeros(Nmesh)
for t in range(Nmesh):
    dy[t] = y[t][0] * 2.
    dy_str.append('dy = ' + str("{:.2e}".format(dy[t])))
    
# Linear coefficient of the L2 norm
q_analytical = math.log(L2_analytical[-1] / L2_analytical[0]) / math.log(dy[-1] / dy[0])

# Analytical sol. built for plot (continuous)
# -------------------------------------------

# Column names
# 0. y
# 1. eta
# 2. u
# 3. v

eta_exact = np.insert(np.loadtxt(file_exact, delimiter=' ', skiprows=1)[:,1], 0, 0.) # Note that null value is added to take into account wall
u_exact = np.insert(np.loadtxt(file_exact, delimiter=' ', skiprows=1)[:,2], 0, 0.)
v_exact = np.insert(np.loadtxt(file_exact, delimiter=' ', skiprows=1)[:,3], 0, 0.)

# Crop arrays for better plot
eta_max = 15.
eta_crp = []
ustar_crp = []
vstar_crp = []

# Numerical sol.
for t in range(Nmesh):
    for i in range(y[t].size):
        if eta[t][i] > eta_max:
            p = int(i - 1)
            break
    eta_crp.append(np.array(eta[t][0:p]))
    ustar_crp.append(np.array(ustar[t][0:p]))
    vstar_crp.append(np.array(vstar[t][0:p]))

# Analytical sol.
for i in range(eta_exact.size):
    if eta_exact[i] > eta_max:
        p = int(i)
        break

eta_exact_crp = np.array(eta_exact[0:p])
u_exact_crp = np.array(u_exact[0:p])
v_exact_crp = np.array(v_exact[0:p])

#-------------------------------------------------------------------------------------------------------------------------------
#                                                         PLOT RESULTS                                                          
#-------------------------------------------------------------------------------------------------------------------------------
# Plot parameters
# ---------------

SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 16
rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the figure title
rc('axes', labelsize=BIGGER_SIZE)     # fontsize of the x and y labels
rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
rc('text', usetex=True)
rc('font', family='serif')
style_plt = ['r-d', 'g-X', 'b-o', 'm-v'] # style use for curve plot

# Figure 1
# --------
plt.clf()
plt.title('Blasius convergence')
plt.plot(dy, L2_analytical, 'b-d', label='Simulation')

#x = np.array([dy[-2], dy[0]])
#y = np.array([1. / dy[-2], 1. / dy[0]])
#plt.plot(x, y, 'k-', label=r'$\propto 1/N^2$')

plt.text(1.3e-4, 1.2e-2, r'$q$={:.3f}'.format(q_analytical), fontsize=MEDIUM_SIZE, color='b')
plt.xlabel(r'$\Delta y_l$')
plt.ylabel(r'$\mathcal{L}^2$')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(b=True, which='major', color='k', linestyle='-')
plt.grid(b=True, which='minor', color='k', linestyle=':')
# ax = plt.gca()
# ax.set_xticks([20, 40, 80], minor=True)
# plt.xticks([20, 40, 80], [20, 40, 80])
plt.savefig('Blasius_L2.pdf', bbox_inches='tight')
# plt.show()

# Figure 2
# --------
plt.clf()
plt.title('Velocity profile x')
for t in range(Nmesh):
    plt.plot(eta_crp[t], ustar_crp[t], style_plt[t], label=dy_str[t])
plt.plot(eta_exact_crp, u_exact_crp, 'k-', label='Analytical')
plt.xlabel(r'$\eta$')
plt.ylabel(r'$u^{*}$')
plt.legend()
plt.savefig('Blasius_vx.pdf', bbox_inches='tight')

# Figure 3
# --------
plt.clf()
plt.title('Velocity profile y')
for t in range(Nmesh):
    plt.plot(eta_crp[t], vstar_crp[t], style_plt[t], label=dy_str[t])
plt.plot(eta_exact_crp, v_exact_crp, 'k-', label='Analytical')
plt.xlabel(r'$\eta$')
plt.ylabel(r'$v^{*}$')
plt.legend()
plt.savefig('Blasius_vy.pdf')