#!/usr/bin/python3

#-------------------------------------------------------------------------------------------------------------------------------
#                                                           OBJECTIVE                                                           
#-------------------------------------------------------------------------------------------------------------------------------
# This script is for the post-treatment of a supersonic viscous flow over a heated flat plate.
# The interaction with a flat plate at fixed temperature produces thermal and dynamic boundary layers
# as well as a Leading Edge Interaction Shock (LEIS).
# Analytical solutions are not avalaible, thus only a comparison with numerical results is done here.
# Comparison is done with results of following publications:
# - Thevand, N., Daniel, E., & Loraud, J. C. (1999). On high‐resolution schemes for solving unsteady compressible two‐phase 
#   dilute viscous flows. International journal for numerical methods in fluids, 31(4), 681-702.
# - Jacobs, P. A. (1991). Single-block Navier-Stokes integrator.
#   INSTITUTE FOR COMPUTER APPLICATIONS IN SCIENCE AND ENGINEERING HAMPTON VA.

#-------------------------------------------------------------------------------------------------------------------------------
#                                                            IMPORTS                                                            
#-------------------------------------------------------------------------------------------------------------------------------
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------------------------------------------------------
#                                                        PARAMETERS TO DEFINE                                           
#-------------------------------------------------------------------------------------------------------------------------------
# Observables
# -----------

# Single file extracted from ParaView
f_ecogen = 'ecogen.csv'
# 0. T
# 1. ux
# 2. uy
# 3. uz
# 4. x
# 5. y
# 6. z

# Two files from digitization
f_ux_jacobs = 'ux_jacobs.txt'
# 0. y
# 1. ux
f_temp_jacobs = 'temp_jacobs.txt'
# 0. y
# 1. T

# Single file from Daniel's code output
f_daniel = 'daniel.txt'
# 0. y
# 1. ux
# 2. uy
# 3. p
# 4. T
# 5. rho
# Other columns are not useful

# Parameter
ymax = 0.02 # Max height for plot

#-------------------------------------------------------------------------------------------------------------------------------
#                                                        READING RESULTS                                               
#-------------------------------------------------------------------------------------------------------------------------------
# Reading input files
# -------------------
y_ecogen = np.loadtxt(f_ecogen, delimiter=',', skiprows=1)[:,5]
ux_ecogen = np.loadtxt(f_ecogen, delimiter=',', skiprows=1)[:,1]
temp_ecogen = np.loadtxt(f_ecogen, delimiter=',', skiprows=1)[:,0]

y_ux_jacobs = np.loadtxt(f_ux_jacobs, delimiter=',', skiprows=1)[:,0]
ux_jacobs = np.loadtxt(f_ux_jacobs, delimiter=',', skiprows=1)[:,1]

y_temp_jacobs = np.loadtxt(f_temp_jacobs, delimiter=',', skiprows=1)[:,0]
temp_jacobs = np.loadtxt(f_temp_jacobs, delimiter=',', skiprows=1)[:,1]

y_daniel = np.loadtxt(f_daniel, delimiter=',')[:,0]
ux_daniel = np.loadtxt(f_daniel, delimiter=',')[:,1]
temp_daniel = np.loadtxt(f_daniel, delimiter=',')[:,4]

# Crop results
# ------------
index = 0

for i in range(y_ecogen.size):
    if y_ecogen[i] > ymax:
        index = int(i)
        y_ecogen = np.array(y_ecogen[0:index])
        ux_ecogen = np.array(ux_ecogen[0:index])
        temp_ecogen = np.array(temp_ecogen[0:index])
        break

for i in range(y_ux_jacobs.size):
    if y_ux_jacobs[i] > ymax:
        index = int(i)
        y_ux_jacobs = np.array(y_ux_jacobs[0:index])
        ux_jacobs = np.array(ux_jacobs[0:index])
        break

for i in range(y_temp_jacobs.size):
    if y_temp_jacobs[i] > ymax:
        index = int(i)
        y_temp_jacobs = np.array(y_temp_jacobs[0:index])
        temp_jacobs = np.array(temp_jacobs[0:index])
        break

for i in range(y_daniel.size):
    if y_daniel[i] > ymax:
        index = int(i)
        y_daniel = np.array(y_daniel[0:index])
        ux_daniel = np.array(ux_daniel[0:index])
        temp_daniel = np.array(temp_daniel[0:index])
        break

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

# Figure Velocity on x
# --------------------
plt.clf()
plt.title('Velocity profile')

plt.plot(y_ecogen, ux_ecogen, 'k-d', label='ECOGEN')
plt.plot(y_ux_jacobs, ux_jacobs, 'b-o', label='Jacobs')
plt.plot(y_daniel, ux_daniel, 'r-x', label='Daniel')

plt.xlabel(r'$y \ (m)$')
plt.ylabel(r'$u \ (m.s^{-1})$')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('velocity.pdf', pad_inches=0.02)

# Figure temperature
# ------------------
plt.clf()
plt.title('Temperature profile')

plt.plot(y_ecogen, temp_ecogen, 'k-d', label='ECOGEN')
plt.plot(y_temp_jacobs, temp_jacobs, 'b-o', label='Jacobs')
plt.plot(y_daniel, temp_daniel, 'r-x', label='Daniel')

plt.xlabel(r'$y \ (m)$')
plt.ylabel(r'$T \ (K)$')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('temperature.pdf', pad_inches=0.02)
