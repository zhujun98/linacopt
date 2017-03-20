#!/usr/bin/python
"""
This is a basic example showing how to optimize the emittance in ASTRA
with a local search optimizer.

The solution is 0.141 at laser_spot = 0.040 and main_sole_b = 0.2815.
"""
import sys, os
sys.path.append(os.path.expanduser('~') + "/myscripts/linac_opt/")

from linac_opt import LinacOpt

# ----------------------------
# Instantiate the optimization

opt_test = LinacOpt(path_name='./astra_basic',
                    input_file='injector.in',
                    input_template='injector.in.000',
                    particle_type='astra')

# --------------------
# Set up the optimizer
opt_test.set_optimizer('sdpen')
opt_test.optimizer.setOption('alfa_stop', 1e-2)

# --------------
# Add fit points
opt_test.fit_points.set_point('out', 'injector.0400.001')

# -------------
# Add objective
def f1(fits):
    print "Horizontal emittance: {:.4f} um".format(fits.out.emitx*1.e6)
    return fits.out.emitx*1.e6
opt_test.opt_prob.set_obj('emitx_um', f1)

# ---------------
# Add constraints
def g1(fits):
    return fits.out.n0
opt_test.opt_prob.set_con('npar', g1, equal=500)

# ------------------------------------------------
# Add variables, co-variables and static-variables
opt_test.opt_prob.set_var('laser_spot', value=0.1, lower=0.04, upper=0.3)
opt_test.opt_prob.set_var('main_sole_b', value=0.1, lower=0.0, upper=0.4)

# --------------------
# Run the optimization
opt_test.solve('astra')
# opt_test.solve('mpirun -np 2 astra_r62_Linux_x86_64_OpenMPI_1.6.1')
