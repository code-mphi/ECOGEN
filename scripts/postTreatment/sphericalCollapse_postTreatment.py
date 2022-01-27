#-------------------------------------------------------------------------------------------------------------------------------
#                                                           OBJECTIVE                                                           
#-------------------------------------------------------------------------------------------------------------------------------
# This script is for the post-treatment of spherical bubble collapse.
# Indeed, it allows to compute the radius over time of the bubble and compare it to the semi-analytical solution of the
# compressible Keller--Miksis equation.
# To do so, it is necessary to define the corresponsding parameters and data to post-treat in the following section. The data
# must contain the total volume of the bubble in the "infoCalcul.out" file, which is not a feature always exported. Appropriate
# commented sections in Run::initialize(), Run::solver() and Output::saveInfos() methods must therefore be uncommented. The
# concerned variable is m_volumePhaseK.

#-------------------------------------------------------------------------------------------------------------------------------
#                                                            IMPORTS                                                            
#-------------------------------------------------------------------------------------------------------------------------------
import math
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

#-------------------------------------------------------------------------------------------------------------------------------
#                                           PARAMETERS TO DEFINE / INITIAL CONDITIONS                                           
#-------------------------------------------------------------------------------------------------------------------------------
# Dimensionful variables
# ----------------------
rhoL = 1000.
cL = 1631.652
gamma = 1.4
muL = 0.
sigma = 0.
Pv = 0. #Pv = 1500.; #Pv = 3550.
Pinf = 50.6625e5
R0 = 1.e-4
Pb0 = 3550. + 2.*sigma/R0
P = 0.
w = 0.
k = (Pb0 - Pv) * R0**(3.*gamma)

# Number of points for the semi-analytical Keller--Miksis solution
# ----------------------------------------------------------------
time_pts = 2000

# Load data in files
# ------------------
file_1D = open("../../results/pressureVelocityEq1DsphericalCollapse_Pratio1427_100cellsPerD/infoCalcul.out", "r")
file_2D = open("../../results/pressureVelocityEq2DsphericalCollapse_Pratio1427_200cellsPerD/infoCalcul.out", "r")

#-------------------------------------------------------------------------------------------------------------------------------
#                                                         COMPUTATIONS                                                          
#-------------------------------------------------------------------------------------------------------------------------------
# Dimensionless variables
# -----------------------
x_init = R0/R0
P0 = Pinf - Pv
k_dimLess = k / (P0 * R0**(3.*gamma))
sigma_dimLess = 2.*sigma/R0/P0
cL_dimLess = cL*(rhoL/P0)**(1./2.)
nuL = 4.*muL/R0/(rhoL*P0)**(1./2.)
w_dimLess = w*R0*(rhoL/P0)**(1./2.)

# Radius and interface velocity
# -----------------------------
u_star_eq = 0.
u_star_outEq = -(1.-Pb0/Pinf)/(cL/(Pinf/rhoL)**0.5)
x0_eq = [x_init, u_star_eq]
x0_outEq = [x_init, u_star_outEq]

# Times
# -----
t_TC = 0.915*(rhoL/(Pinf-Pv))**(1./2.)*R0
t_dimLessCoeff = (P0/rhoL)**(1./2.)/R0
t0 = 0.
t_end = 2.*t_TC*t_dimLessCoeff
tspan = np.linspace(t0, t_end, time_pts) #The points of evaluation of solution
tspan_plot = tspan/t_dimLessCoeff/t_TC

# Evaluate the Keller-Miksis ODE with substitution parameter to obtain a set of first order equations
# ---------------------------------------------------------------------------------------------------
def MK(t_dimLess, x, P0, k_dimLess, gamma, nuL, sigma_dimLess, cL_dimLess, P, w_dimLess):
    dxdt = np.zeros(2)
    dxdt[0] = x[1]
    dxdt[1] =  ( (dxdt[0]**3.)/2. + dxdt[0]*((1.-3.*gamma)*k_dimLess*x[0]**(-3.*gamma) - 1.)
                - cL_dimLess*(3./2.*dxdt[0]**2. + nuL*dxdt[0]/x[0] + sigma_dimLess/x[0] - k_dimLess*x[0]**(-3.*gamma) + 1.)
                + (1. + dxdt[0]/cL_dimLess)*P*cL_dimLess/P0*math.sin(w_dimLess*(t_dimLess + x[0]/cL_dimLess)) ) / (nuL - x[0]*(dxdt[0] - cL_dimLess))
    return dxdt

# Solve the compressible Keller-Miksis equation
# ---------------------------------------------
# Equilibrium
x_eq = np.zeros((len(tspan), len(x0_eq)))   #Array for solution (R, u)
x_eq[0, :] = x0_eq
r_eq = integrate.ode(MK).set_integrator("dopri5")  #Choice of method
r_eq.set_initial_value(x0_eq, t0).set_f_params(P0, k_dimLess, gamma, nuL, sigma_dimLess, cL_dimLess, P, w_dimLess)   #Initial values (y, t) and then *args
for i in range(1, tspan.size):
    x_eq[i, :] = r_eq.integrate(tspan[i]) #Get one more value, add it to the array
    if not r_eq.successful():
        raise RuntimeError("Could not integrate")
# Out of equilibrium
x_outEq = np.zeros((len(tspan), len(x0_outEq)))   #Array for solution (R, u)
x_outEq[0, :] = x0_outEq
r_outEq = integrate.ode(MK).set_integrator("dopri5")  #Choice of method
r_outEq.set_initial_value(x0_outEq, t0).set_f_params(P0, k_dimLess, gamma, nuL, sigma_dimLess, cL_dimLess, P, w_dimLess)   #Initial values (y, t) and then *args
for i in range(1, tspan.size):
    x_outEq[i, :] = r_outEq.integrate(tspan[i]) #Get one more value, add it to the array
    if not r_outEq.successful():
        raise RuntimeError("Could not integrate")

#-------------------------------------------------------------------------------------------------------------------------------
#                                               SIMULATION RESULTS DETERMINATION                                                
#-------------------------------------------------------------------------------------------------------------------------------
# Radius computation
# ------------------
# 1D - volume
# 2D - (volume*4/math.pi)**(1/2)     ; note that we only compute 1/4th of the disk
# 3D - (volume*8*3/4/math.pi)**(1/3) ; note that we only compute 1/8th of the sphere

# 1D
# --
lines_1D = file_1D.readlines()
del lines_1D[0] #Delete the first line of the file (useless)
time_1D = []
volume_1D = []
for line in lines_1D:
    split_line = line.split()
    volume_1D.append(float(split_line[7])) #The volume is the 8th variable in the infoCalcul.out file
    time_1D.append(float(split_line[2]))   #The time is the 3rd variable in the infoCalcul.out file
R0_1D = volume_1D[0]
non_dim_radius_1D = [element /R0_1D for element in volume_1D]     #Non-dimensionalized radius
non_dim_time_1D = [element /t_TC*R0/R0_1D for element in time_1D] #Non-dimensionalized time

# 2D
# --
lines_2D = file_2D.readlines()
del lines_2D[0] #Delete the first line of the file (useless)
time_2D = []
volume_2D = []
for line in lines_2D:
    split_line = line.split()
    volume_2D.append(float(split_line[7])) #The volume is the 8th variable in the infoCalcul.out file
    time_2D.append(float(split_line[2]))   #The time is the 3rd variable in the infoCalcul.out file
R0_2D = (volume_2D[0]*4/math.pi)**(1/2)
non_dim_radius_2D = [(element *4/math.pi)**(1/2) / R0_2D for element in volume_2D] #Non-dimensionalized radius
non_dim_time_2D = [element /t_TC*R0/R0_2D for element in time_2D]                  #Non-dimensionalized time

# 3D
# --
# lines_3D = file_3D.readlines()
# del lines_3D[0] #Delete the first line of the file (useless)
# time_3D = []
# volume_3D = []
# for line in lines_3D:
#     split_line = line.split()
#     volume_3D.append(float(split_line[7])) #The volume is the 8th variable in the infoCalcul.out file
#     time_3D.append(float(split_line[2]))   #The time is the 3rd variable in the infoCalcul.out file
# R0_3D = (volume_3D[0]*8*3/4/math.pi)**(1/3)
# non_dim_radius_3D = [(element *8*3/4/math.pi)**(1/3) / R0_3D for element in volume_3D] #Non-dimensionalized radius
# non_dim_time_3D = [element /t_TC*R0/R0_3D for element in time_3D]                      #Non-dimensionalized time

#-------------------------------------------------------------------------------------------------------------------------------
#                                                         PLOT RESULTS                                                          
#-------------------------------------------------------------------------------------------------------------------------------
# Plot parameters
# ---------------
SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 16
rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the figure title
rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
rc('text', usetex=True)
rc('font', family='serif')

# Figure 1
# --------
plt.title('Axi-symmetry validation')
plt.plot(tspan_plot, x_outEq[:, 0], 'k-', label='Keller--Miksis, $\dot{R}_0 \\neq 0$')
# plt.plot(tspan_plot, x_eq[:, 0], 'k:', label='Keller--Miksis, $\dot{R}_0 = 0$')
plt.plot(non_dim_time_1D, non_dim_radius_1D, 'r--', label='1D - 100c/d')
plt.plot(non_dim_time_2D, non_dim_radius_2D, 'b-.', label='2D - 200c/d')
plt.xlabel(r'$t/t_c$')
plt.ylabel(r'$R/R_0$')
plt.legend()
# plt.grid(True)
plt.savefig('fig1.pdf')
plt.show()
