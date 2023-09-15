#!/usr/bin/env python3

#-------------------------------------------------------
#                       OBJECTIVE                       
#-------------------------------------------------------

# This script allows to compute the saturation state using a given
# pressure or temperature for a liquid/gas mixture.
# The Stiffened Gas or Ideal Gas (pinf = 0) equation of state (EOS) 
# is used for each phase.

#-------------------------------------------------------
#                       IMPORTS
#-------------------------------------------------------
from math import exp, log

#-------------------------------------------------------
#                 PARAMETERS TO DEFINE
#-------------------------------------------------------

# EOS parameters
# Liquid
cvl = 1816.    # Heat capacity at constant volume   (J.kg-1.K-1)
gammal = 2.35  # Adiabatic index                    (-)
pinfl = 1.e9   # Matter cohesion parameter          (Pa)
ql = -1167.e3  # Energy constant                    (J.kg-1)
qpl = 0.       # Entropy constant                   (J.kg-1.K-1)

# Gas
cvg = 1040.    # Heat capacity at constant volume   (J.kg-1.K-1)
gammag = 1.43  # Adiabatic index                    (-)
pinfg = 0.     # Matter cohesion parameter          (Pa)
qg = 2030.e3   # Energy constant                    (J.kg-1)
qpg = -23.4e3  # Entropy constant                   (J.kg-1.K-1)

# Input state condition
p = 100.e5     # Saturation pressure    (Pa)
T = 373.       # Saturation temperature (K)
computeSatFromTemp = False # Choice on the method to compute saturation state

#-------------------------------------------------------
#                     COMPUTATIONS
#-------------------------------------------------------

cpl = gammal * cvl # Heat capacity at constant pressure (J.kg-1.K-1)
cpg = gammag * cvg # Heat capacity at constant pressure (J.kg-1.K-1)

# Compute Psat coefficients 
A = (cpl - cpg + qpg - qpl) / (cpg - cvg)
B = (ql - qg) / (cpg - cvg)
C = (cpg - cpl) / (cpg - cvg)
D = (cpl - cvl) / (cpg - cvg)

print('Running...')

err = 1.
it = 0

if (computeSatFromTemp):
  print('Compute saturation state from temperature')
  psat1 = 2.e5 # Initial guess
  while (abs(err) > 1.e-5):
    f = psat1 + pinfg - exp( A + B / T + C * log(T) ) * (psat1 + pinfl)**D
    df = 1. - exp( A + B / T + C * log(T)) * D * (psat1 + pinfl)**(D - 1.)

    psat2 = psat1 - f / df
    err = (psat2 - psat1) / (0.5 * (psat1 + psat2))
    print('Iteration:', it, '| Psat =', psat2, '| Psat(old) =', psat1)
    psat1 = psat2
    
    it += 1
    if (it > 50):
      print('WARNING: Newton-Raphson has not converged.')
  # Save final pressure
  p = psat1
  print('Temp =', T, 'K')
  print('Psat =', p, 'Pa')

else:
  print('Compute saturation state from pressure')
  Tsat1 = 350. # Initial guess
  while (abs(err) > 1.e-5):
    f = log(p + pinfg) - A - B / Tsat1 - C * log(Tsat1) - D * log(p + pinfl)
    df =  B / (Tsat1**2.) - C / Tsat1

    Tsat2 = Tsat1 - f / df
    err = (Tsat2 - Tsat1) / (0.5 * (Tsat1 + Tsat2))
    print('Iteration:', it, '| Tsat =', Tsat2, '| Tsat(old) =', Tsat1)
    Tsat1 = Tsat2
    
    it += 1
    if (it > 50):
      print('WARNING: Newton-Raphson has not converged.')
  # Save final temperature
  T = Tsat1
  print('p =', p, 'Pa')
  print('Tsat =', T, 'K')

# Compute saturation phase densities
rhog = (p+pinfg) / (gammag - 1.) / cvg / T
rhol = (p+pinfl) / (gammal - 1.) / cvl / T

#-------------------------------------------------------
#                   WRITE RESULTS
#-------------------------------------------------------

print('rhoGSat =', rhog, 'kg/m3')
print('rhoLSat =', rhol, 'kg/m3')

print('Done.')