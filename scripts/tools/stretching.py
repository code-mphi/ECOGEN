# ----------------------
#       OBJECTIVE
# ----------------------

# Compute the number of mesh cells to use when stretching along a direction
# and given the length of the domain, the geometric ratio and the length
# of the first cell.
# Also compute the number of mesh cells to use for a refined or coarse mesh
# where the length of the first cell is half or twice of the one of the
# previous mesh.

# Example:
# A 1D mesh has Nx cells with a geometric progression 'q' along x-direction.
# The length of first cell is dx0.
# One would like to run a simulation with half length of first cell (dx0/2)
# using the same geometric progression. How to choose the new Nx?

# Recall:
#  Geometric sequence is defined by:
#   u_n+1 = u_n * q
#   u_n = u_0 * q^n
#  The sum is given by: sum_{k = 0}^{n} u_k = u_0 (1 - q^n) / (1 - q)

# ----------------------
#       IMPORTS
# ----------------------
import math

# ----------------------
#       PARAMETERS
# ----------------------
h = 0.3    # Length of the domain in the stretching direction
q = 1.045  # Geometric ratio
dx0_initialMesh = 6.9e-05               # Length of 1st cell of initial mesh
dx0_refinedMesh = 0.5 * dx0_initialMesh # Length of 1st cell of refined mesh
dx0_coarseMesh = 2. * dx0_initialMesh   # Length of 1st cell of coarse mesh

# ----------------------
#       COMPUTATION
# ----------------------

# Geometric progression
# ---------------------
#  Initial mesh
#  ------------
#  x_Ni = x0i * q^Ni
#  Ni: number of cells of initial mesh
#  x0i: length of 1st cell of initial mesh

Ni = round(math.log( 1. - (h / dx0_initialMesh + 1.) * (1. - q) ) / math.log(q) )

# Refined mesh
# ------------
#  x_Nr = x0r * q^Nr
#       = 0.5 * x0i * q^Nr
#  Nr: number of cells of refined mesh (to find)
#  x0r: length of 1st cell of refined mesh

# Nr = round(math.log( 1. -  h / dx0_refinedMesh       * (1. - q) ) / math.log(q) )
Nr = round(math.log( 1. - (h / dx0_refinedMesh + 1.) * (1. - q) ) / math.log(q) )

# Coarse mesh
# -----------
#  x_Nc = x0c * q^Nc
#       = 2. * x0i * q^Nc
#  Nc: number of cells of coarse mesh (to find)
#  x0c: length of 1st cell of coarse mesh

Nc = round(math.log( 1. - (h / dx0_coarseMesh + 1.) * (1. - q) ) / math.log(q) )

# ----------------------
#         OUTPUT
# ----------------------
print('------------ Initial mesh ------------')
print('First cell length dx = {}'.format(dx0_initialMesh))
print('Number of cells to use N = {}'.format(Ni))
print('------------ Refined mesh ------------')
print('New first cell length dx = {}'.format(dx0_refinedMesh))
print('Number of cells to use N = {}'.format(Nr))
print('------------ Coarse mesh -------------')
print('New first cell length dx = {}'.format(dx0_coarseMesh))
print('Number of cells to use N = {}'.format(Nc))