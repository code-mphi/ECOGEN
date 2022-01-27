#!/usr/bin/python3

#-------------------------------------------------------------------------------------------------------------------------------
#                                                           OBJECTIVE                                                           
#-------------------------------------------------------------------------------------------------------------------------------
# This script (1/2) is for the post-treatment of a subsonic flow over a flat flate, i.e. Blasius' test case.
# It aims to compare steady dimensionless normal/tangential velocities with the analytical solution.
# The analytical solution is based on the book:
#  'I do like CFD', vol. 1, Katate Masatsuka
# The comparison is done through a convergence study with L2 norms.

#-------------------------------------------------------------------------------------------------------------------------------
#                                                            IMPORTS                                                            
#-------------------------------------------------------------------------------------------------------------------------------
import math
import numpy as np

#-------------------------------------------------------------------------------------------------------------------------------
#                                           PARAMETERS TO DEFINE / INITIAL CONDITIONS                                           
#-------------------------------------------------------------------------------------------------------------------------------
# Variables
# ---------
uinf = 102.    # Free stream velocity
rho = 1.2      # Density (supposed constant for Blasius solution)
mu = 1.82e-5   # Dynamic viscosity
xplate = 0.02  # Position of leading edge plate
x = 0.405      # Horizontal position for data extraction
h = 0.01       # Max. vertical position

# Observables
# -----------
files = ['mesh_200x90.csv', 'mesh_400x105.csv', 'mesh_800x120.csv', 'mesh_1600x135.csv'] # Filenames (under current directory)

# Parameters for reference analytical solution
# --------------------------------------------
dh0 = 1.e-4   # First vertical position for plot analytical sol.
q = 1.045     # Geometric sequence ratio for vertical position

#-------------------------------------------------------------------------------------------------------------------------------
#                                                         COMPUTATIONS                                                          
#-------------------------------------------------------------------------------------------------------------------------------
nu = mu/rho             # Kinematic viscosity
xref = x - xplate       # Reference position from the beginning of leading edge
Re = uinf * xref / nu   # External Reynolds number

Ny = round(math.log(h / dh0) / math.log(q)) # Number of vertical position for plot exact sol.

#-------------------------------------------------------------------------------------------------------------------------------
#                                                            FUNCTIONS                                                            
#-------------------------------------------------------------------------------------------------------------------------------
def rhs(F):
    """Build right hand side of the ODE eq. 7.15.53.

    This function is used during ODE integration insde blasius() function.

    Args:
        F: vector F to integrate (3 components)

    Returns:
        G: RHS vector G (3 components)
    """
    
    f = F[0]
    f1 = F[1]
    f2 = F[2]

    G = np.array([f1, f2, - 0.5 * f * f2])
    return G

def blasius(etap, Re):
    """Compute the Blasius analytical solution (u,v) at a given position (xp, yp).
    
    This method uses the one-step integration given in the book 'I do like CFD',
    vol. 1, Katate Masatsuka (see section 'Flat Plate Boundary Layer').

    Args:
        etap: dimensionless Blasius variable
        Re:   free stream Reynolds number at xref position

    Returns:
        u: dimensionless horizontal velocity
        v: dimensionless vertical velocity
    """
    
    F = np.zeros(3)           # Left hand side of ODE (see eq. 7.15.53)
    K1 = np.zeros(3)          # Runge--Kutta variable K1 (3 components vector)
    K2 = np.zeros(3)          # Runge--Kutta variable K2
    K3 = np.zeros(3)          # Runge--Kutta variable K3
    
    f2_0 = 0.3320573362151946 # Pre-computed initial value (see eq. 7.15.59)
    deta = 1.e-3              # Increment for ODE integration
    done = False
    one_over_three = 1. / 3.
    one_over_six = 1. / 6.

    # Integrate ODE to eta = etap by RK4 scheme
    # 1. Initial values
    eta = 0.
    F[0] = 0.
    F[1] = 0.
    F[2] = f2_0
    
    # 2. Steps to eta = etap
    if etap > 0.:
        while True:
            if (eta + deta > etap):
                deta = etap - eta
                done = True
            eta += deta
            K1 = F + 0.5 * deta * rhs(F)
            K2 = F + 0.5 * deta * rhs(K1)
            K3 = F + deta * rhs(K2)
            F = (K1 + 2. * K2 + K3 - F) * one_over_three + deta * rhs(K3) * one_over_six

            if (done):
                break

    # Solution at eta = etap, i.e. (x,y) = (xp,yp)
    f = F[0]
    f1 = F[1]
    u = f1 
    v = 0.5 / math.sqrt(Re) * (etap * f1 - f)
    return u, v

#-------------------------------------------------------------------------------------------------------------------------------
#                                               SIMULATION RESULTS DETERMINATION                                                
#-------------------------------------------------------------------------------------------------------------------------------
# Results
# -------
Nmesh = len(files)  # Number of test cases
y = []              # Array of vertical positions
u = []              # Array of normal velocities
v = []              # Array of tangential velocities
ustar = []          # Array of dimensionless normal velocity
vstar = []          # Array of dimensionless tangential velocity
eta = []            # Array of Blasius dimensionless variable
u_bls = []          # Array of dimensionless analytical normal velocity
v_bls = []          # Array of dimensionless analytical tangential velocity

# Read data of each mesh
# Column names
# 0. Vx
# 1. Vy
# 2. Vz
# 3. x
# 4. y
# 5. z
    
count = 0
for f in files:
    y.append(np.loadtxt(f, delimiter=',', skiprows=1)[:,4])
    u.append(np.loadtxt(f, delimiter=',', skiprows=1)[:,0])
    v.append(np.loadtxt(f, delimiter=',', skiprows=1)[:,1])
    
    ustar.append(np.zeros(y[count].size))
    vstar.append(np.zeros(y[count].size))
    eta.append(np.zeros(y[count].size))
    u_bls.append(np.zeros(y[count].size))
    v_bls.append(np.zeros(y[count].size))    
    count += 1

# L2 convergence with ponderation by cells' height
dy = []
for t in range(Nmesh):
    dy.append(np.zeros(y[t].size))
    dy[t][0] = 2. * y[t][0]
    
for t in range(Nmesh):
    sumDy = 0.
    for i in range(1, y[t].size):
        sumDy += dy[t][i-1] # Current state
        dy[t][i] = 2. * (y[t][i] - sumDy)


for t in range(Nmesh):
    for i in range(y[t].size):
        eta[t][i] = y[t][i] / xref * math.sqrt(Re)
        ustar[t][i] = u[t][i] / uinf
        vstar[t][i] = v[t][i] / uinf # * math.sqrt(Re)    
    
    # Analytical solutions for L2 norm
    # --------------------------------
    print('Computing analytical solution for L2 norm of file {}...'.format(files[t]))
    for i in range(y[t].size):
        print(' {}%'.format(int(i / y[t].size * 100)))
        u_bls[t][i], v_bls[t][i] = blasius(eta[t][i], Re)

# Convergence
# -----------
# With analytical solution being surrogate truth
L2_analytical = np.zeros(Nmesh) # L2 norm
for t in range(Nmesh):
    num = 0.
    den = 0.
    for i in range(y[t].size):
        # num += (u_bls[t][i] - ustar[t][i]) ** 2.
        num += dy[t][i] * (u_bls[t][i] - ustar[t][i]) ** 2.        
        den += dy[t][i] 
    # L2_analytical[t] = (num / y[t].size) ** 0.5
    L2_analytical[t] = (num / den) ** 0.5
        
# Analytical solution for plot
# ----------------------------
u_plt = np.zeros(Ny)    # Array of dimensionless analytical normal velocity
v_plt = np.zeros(Ny)    # Array of dimensionless analytical tangential velocity
y_plt = np.fromiter((dh0 * q**i for i in range(Ny)), float) # Array of vertical position (following geometric sequence)
eta_plt = np.fromiter((yi / xref * math.sqrt(Re) for yi in y_plt), float) # Array of dimensionless position, i.e. Blasius var.

print('Computing analytical solution for plot...')
for i in range(Ny):
    print(' {}%'.format(int(i / Ny * 100)))
    eta_plt[i] = y_plt[i] / xref * math.sqrt(Re)
    u_plt[i], v_plt[i] = blasius(eta_plt[i], Re)

#-------------------------------------------------------------------------------------------------------------------------------
#                                                         WRITE RESULTS                                                          
#-------------------------------------------------------------------------------------------------------------------------------

# Results for simulations
# -----------------------
print('Writing results...')
# Each simulation has its own output file
columns = 'y' + ' ' + 'eta' + ' ' + 'ustar' + ' ' + 'vstar' + ' ' \
    + 'u_bls' + ' ' + 'v_bls' + ' ' + 'L2_analytical' + '\n'

files = [files[t][:-4] for t in range(Nmesh)] # remove csv extension

for t in range(Nmesh):
    file = 'res_' + files[t] + '.txt'
    print('Writing file {} done'.format(file))
    data = np.array([y[t], eta[t], ustar[t], vstar[t], 
                     u_bls[t], v_bls[t], L2_analytical[t]], dtype=object)
    data = data.T
    with open(file, 'w') as f:
        f.write(columns)
        for i in range(y[t].size):
            f.write(str(y[t][i]) + ' '  
                    + str(eta[t][i]) + ' ' 
                    + str(ustar[t][i]) + ' ' 
                    + str(vstar[t][i]) + ' '
                    + str(u_bls[t][i]) + ' '
                    + str(v_bls[t][i]) + ' '
                    + str(L2_analytical[t])+ '\n')
            
# Results for plotting analytical sol.
# ------------------------------------
file = 'exact.txt'
print('Writing file {} done'.format(file))
columns = 'y' + ' ' + 'eta' + ' ' + 'u' + ' ' + 'v' + '\n'
with open(file, 'w') as f:
    f.write(columns)
    for i in range(Ny):
        f.write(str(y_plt[i]) + ' '
                + str(eta_plt[i]) + ' '
                + str(u_plt[i]) + ' '
                + str(v_plt[i]) + '\n')

print('')        
print('Results are stored in files mentionned above.\nTo plot use script blasius_plots.py.')