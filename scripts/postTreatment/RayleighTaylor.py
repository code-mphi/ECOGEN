#-------------------------------------------------------------------------------------------------------------------------------
#                                                           OBJECTIVE                                                           
#-------------------------------------------------------------------------------------------------------------------------------
# This script is for the post-treatment of the Rayleigh--Taylor instability.
# It aims to compare the interface front displacement of the simulations with the analytical inviscid and/or viscous solutions.
# The analytical solutions are based on the papers of:
#  Nourgaliev, R. R., Liou, M. S., & Theofanous, T. G. (2008). Numerical prediction of interfacial instabilities: Sharp interface method (SIM). Journal of Computational Physics, 227(8), 3940-3970.
#  Duff, R. E., Harlow, F. H., Hirt, C. W. (1962). Effects of diffusion on interface instability between gases. The Physics of Fluids, 5(4), 417-425.
#  Cook, A. W., Cabot, W., & Miller, P. L. (2004). The mixing transition in Rayleigh-Taylor instability. Journal of Fluid Mechanics, 511, 333.
# The comparison is done through a convergence study with L2 norms.

#-------------------------------------------------------------------------------------------------------------------------------
#                                                            IMPORTS                                                            
#-------------------------------------------------------------------------------------------------------------------------------
import math
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------------------------------------------------------
#                                           PARAMETERS TO DEFINE / INITIAL CONDITIONS                                           
#-------------------------------------------------------------------------------------------------------------------------------
# Variables
# ---------
rho_l = 1.         #Density light fluid
g = 9.81           #Gravity
wave_length = 0.2  #Wave-length (characteristic length)
At = 0.5           #Atwood number
R = 76.            #Reynolds number
m = 3.             #Viscosity ratio. Note that herein, nu_l = nu_h (h for heavy fluid). Therefore, this ratio is the density ratio

# Observables
# -----------
Ntest = 4                                      #Number of test cases
mesh = np.array([10, 20, 40, 80])              #Resolution in the wave-length direction
mesh = mesh * 2
epsilon_t = np.zeros(Ntest)                    #Interface front displacement
t0 = 0.05                                      #Initial time (when the flow is already started)
t = 1. - t0                                    #Time
front_position_at_t0 = 0.5                     #Initial interface front position
epsilon_0 = 0.5115 - front_position_at_t0      #Initial interface front displacement
epsilon_t[0] = 0.827 - front_position_at_t0
epsilon_t[1] = 0.7902 - front_position_at_t0
epsilon_t[2] = 0.78234 - front_position_at_t0
epsilon_t[3] = 0.78085 - front_position_at_t0

# Analytical-Rayleigh--Taylor, dimensionless growth factor
# --------------------------------------------------------
K_RT_analytical = 0.28

#-------------------------------------------------------------------------------------------------------------------------------
#                                                         COMPUTATIONS                                                          
#-------------------------------------------------------------------------------------------------------------------------------
rho_h = rho_l * (At + 1.) / (1. - At)           #Density heavy fluid
nu = wave_length * ((wave_length * g)**0.5) / R #Kinematic viscosity
mu_l = rho_l * nu                               #Dynamic viscosity light fluid
mu_h = mu_l * m                                 #Dynamic viscosity heavy fluid
tau_RT = (wave_length / g)**0.5                 #Characteristic time
k = 2. * math.pi / wave_length                  #Wave-number
print('Simulation parameters (SI):')
print('lambda =', wave_length)
print('g      =', g)
print('rho_l  =', rho_l)
print('rho_h  =', rho_h)
print('mu_l   =', mu_l)
print('mu_h   =', mu_h)

#-------------------------------------------------------------------------------------------------------------------------------
#                                               SIMULATION RESULTS DETERMINATION                                                
#-------------------------------------------------------------------------------------------------------------------------------
# Results
# -------
t_star = t / tau_RT                #Scaled time
k_RT = k / (g / nu**2.)**(1. / 3.) #Dimensionless wave-number
K_RT = np.zeros(Ntest)             #Dimensionless growth factor
for i in range(0, K_RT.size):
    K_RT[i] = math.log(epsilon_t[i] / epsilon_0) / t / (g**2. / nu)**(1. / 3.)
print('Results:')
print('k_RT =', k_RT)
print('K_RT[-1] =', K_RT[-1])
print('t =', t)
print('Amplitude simu     =', epsilon_t[-1] - epsilon_0)

# Analytical solutions
# --------------------
# Inviscid solution
A_invisc = epsilon_0 * math.exp((At * g * k)**0.5 * t) #Amplitude
print('Amplitude inviscid =', A_invisc)
# Viscous solution
A_viscous = epsilon_0 * math.exp(((At * g * k + nu**2. * k**4.)**0.5 - nu * k**2.) * t) #Amplitude
print('Amplitude viscous  =', A_viscous)
C = 0.85
Ate_over_At = 0.48
alpha = (C / 2.)**2. * Ate_over_At**3.
h = alpha * At * g * t**2. + 2. * (alpha * At * g * epsilon_0)**0.5 * t + epsilon_0
print('alpha =', alpha)
print('h =', h)

# Convergence
# -----------
# With highest resolution being surrogate truth
L2_resolution = np.zeros(Ntest - 1)              #L2 norm
for i in range(0, L2_resolution.size):
    L2_resolution[i] = ((K_RT[i] - K_RT[-1])**2.)**0.5
q_resolution = math.log(L2_resolution[-1] / L2_resolution[0]) / math.log(mesh[0] / mesh[-2]) #Linear coefficient of the L2 norm
# With analytical solution being surrogate truth
L2_analytical = np.zeros(Ntest)                  #L2 norm
for i in range(0, L2_analytical.size):
    L2_analytical[i] = ((K_RT[i] - K_RT_analytical)**2.)**0.5
q_analytical = math.log(L2_analytical[-1] / L2_analytical[0]) / math.log(mesh[0] / mesh[-1]) #Linear coefficient of the L2 norm

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

# Figure 1
# --------
plt.title('Rayleigh--Taylor convergence')
plt.plot(mesh[0:-1], L2_resolution, 'r-d', label='Surrogate truth = highest resolution')
plt.plot(mesh, L2_analytical, 'b-X', label='Surrogate truth = analytical solution')
x = np.array([mesh[0], mesh[2]])
y = np.array([1.6 / mesh[0]**2., 1.6 / mesh[2]**2.])
plt.plot(x, y, 'k-', label=r'$\propto 1/N^2$')
plt.text(40.2, 3.4e-4, r'$q$={:.3f}'.format(q_resolution), fontsize=MEDIUM_SIZE, color='r')
plt.text(40.2, 2.4e-4, r'$q$={:.3f}'.format(q_analytical), fontsize=MEDIUM_SIZE, color='b')
plt.xlabel(r'$N$')
plt.ylabel(r'$\mathcal{L}^2$')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(b=True, which='major', color='k', linestyle='-')
plt.grid(b=True, which='minor', color='k', linestyle=':')
# ax = plt.gca()
# ax.set_xticks([20, 40, 80], minor=True)
# plt.xticks([20, 40, 80], [20, 40, 80])
plt.savefig('RayleighTaylor_L2.pdf')
# plt.show()

