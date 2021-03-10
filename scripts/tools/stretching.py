# ----------------------
#       OBJECTIVE
# ----------------------

# Compute number of mesh cells (using stretching) along a direction 
# to get half or twice length of a 1st cell given mesh

# Example:
# 1D mesh has Nx cells with a geometric progression 'q' along x-direction.
# The length of 1st cell is dx0. 
# One would like to run a simulation with half length of 1st cell (dx0/2) using
# the same geometric progression. How to choose the new Nx? 

# Recall:
# Geometric sequence is defined by: 
#  u_n+1 = u_n * q
#  u_n = u_0 * q^n

#  The sum is given by: sum_{k = 0}^{n} u_k = u_0 (1 - q^n) / (1 - q)

# ----------------------
#       IMPORTS
# ----------------------
import math

# ----------------------
#       PARAMETERS
# ----------------------

# Blasius 
h = 0.3    # Length of the domain in the stretching direction
q = 1.045  # Geometric ratio
N_mesh1 = int(120)
dx0_mesh1 = 6.9e-05         # Length of 1st cell of given mesh
dx0_mesh2 = 0.5 * dx0_mesh1 # Length of 1st cell of wanted mesh

# ----------------------
#       COMPUTATION
# ----------------------

# Geometric progression
# ---------------------
#  Given mesh
#  ----------
#  x_Ng = x0g * q^Ng
#  Ng: number of cells of given mesh
#  x0g: length of 1st cell of given mesh

# Wanted mesh
# -----------
#  x_N2 = x0w * q^Nw 
#       = 0.5 * x0g * q^Nw
#  Nw: number of cells of wanted mesh (to find)
#  x0w: length of 1st cell of wanted mesh

# N = round( math.log( 1. - h / dx0_mesh2 * (1. - q) ) / math.log(q) )
N = round(math.log( 1. - (h / dx0_mesh2 + 1.) * (1. - q) ) / math.log(q) )

# ----------------------
#         OUTPUT
# ----------------------

print('New first cell length dx = {}'.format(dx0_mesh2))
print('Number of cells to use N = {}'.format(N))