#!/usr/bin/env python3

#---------------------------
#         OBJECTIVE         
#---------------------------

# This script allows to compute the numerical saturation
# curve p-v for a given temperature range.
# The liquid and vapor obey SG/IG EOS and are calibrated 
# using the saturation experimental data, for more ref. see
# Le Métayer, O., Massoni, J., & Saurel, R. (2004). 
# Élaboration des lois d'état d'un liquide et de sa vapeur 
# pour les modèles d'écoulements diphasiques. 
# International journal of thermal sciences, 43(3), 265-276.

#---------------------------
#         IMPORTS         
#---------------------------

import numpy as np
from math import floor, exp, log

print('Running...')

#---------------------------
#    PARAMETERS TO DEFINE
#---------------------------

print('Read input EOS parameters...')

# EOS parameters
#  Liquid
cpl = 4183.93       # Heat capacity at constant pressure (J.kg-1.K-1)
cvl = 1479.48       # Heat capacity at constant volume   (J.kg-1.K-1)
pinfl = 8.05163e8   # Matter cohesion parameter          (Pa)
gammal = 2.82798    # Adiabatic index                    (-)
ql = -1.14264e6     # Energy constant                    (J.kg-1)
qpl = 0.            # Entropy constant                   (J.kg-1.K-1)

# Gas
cpg = 1809.74       # Heat capacity at constant pressure (J.kg-1.K-1)
cvg = 1348.71       # Heat capacity at constant volume   (J.kg-1.K-1)
pinfg = 0.          # Matter cohesion parameter          (Pa)
gammag = 1.34183    # Adiabatic index                    (-)
qg = 2.00685e6      # Energy constant                    (J.kg-1)
qpg = -30025.9      # Entropy constant                   (J.kg-1.K-1)

# Temperature range of interest
Tmin = 278.  # Min. temperature (K)
Tmax = 318.  # Max. temperature (K)

#---------------------------
#        COMPUTATIONS
#---------------------------

# Compute Psat coefficients 
A = (cpl - cpg + qpg - qpl) / (cpg - cvg)
B = (ql - qg) / (cpg - cvg)
C = (cpg - cpl) / (cpg - cvg)
D = (cpl - cvl) / (cpg - cvg)

# Temperature range of the calibration (K)
T = np.linspace(Tmin, Tmax, num=50)

psat = np.zeros(T.size)
vl = np.zeros(T.size)
vg = np.zeros(T.size)

print('Compute saturation quantities...')

# Temperature loop
for i in range(T.size):
  progress = int(i / T.size * 100)
  if (progress % 10 == 0):
    print('Progress {}%'.format(progress))
  psat1 = 2.e-5
  err = 1.
  count = 0
  # Compute Psat with Newton-Raphson procedure applied to Gibbs free energy eq.
  while (abs(err) > 1.e-5):
    f = psat1 + pinfg - exp( A + B / T[i] + C * log(T[i]) ) * (psat1 + pinfl)**D
    df = 1. - exp( A + B / T[i] + C * log(T[i])) * D * (psat1 + pinfl)**(D - 1.)
    psat2 = psat1 - f / df
    err = (psat2 - psat1) / (0.5 * (psat1 + psat2))
    psat1 = psat2
    count += 1
    if (count > 50):
      print('WARNING: Newton-Raphson has not converged for temperature ', T[i])
  psat[i] = psat1

  # Compute specific volumes using each phase EOS
  vl[i] = (gammal - 1.) * cvl * T[i] / (psat[i] + pinfl)
  vg[i] = (gammag - 1.) * cvg * T[i] / (psat[i] + pinfg)

#---------------------------
#      WRITE RESULTS
#---------------------------

print('Write results...')

with open('psat-num.out', 'w') as f:
  f.write('T' + ' ' + 'Psat' + ' ' + 'vL' + ' ' + 'vG' + '\n')
  for i in range(T.size):
    f.write(str(T[i]) + ' ' + str(psat[i]) + ' ' + str(vl[i]) + ' ' + str(vg[i]) + '\n')
  f.close()

print('Done.')