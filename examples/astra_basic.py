#!/usr/bin/python
"""
This is a basic example showing how to optimize the emittance in ASTRA
with a local search optimizer.
"""
import sys
import os
sys.path.append(os.path.expanduser('~') + "/myscripts/linac_opt/")

from linac_opt import LinacOpt

# ----------------------------
# Instantiate the optimization

opt_test = LinacOpt(path_name='./astra_basic',
                    input_file='injector.in',
                    input_template='injector.in.000',
                    particle_type='astra',
                    prob_name='opt_test')

# --------------------
# Set up the optimizer
# opt_test.set_optimizer('sdpen')
# opt_test.optimizer.setOption('alfa_stop', 1e-3)

opt_test.set_optimizer('alpso')

opt_test.optimizer.setOption('atol', 1e-2)
opt_test.optimizer.setOption('rtol', 1e-2)
opt_test.optimizer.setOption('SwarmSize', 40)
opt_test.optimizer.setOption('maxOuterIter', 6)
opt_test.optimizer.setOption('maxInnerIter', 3)
opt_test.optimizer.setOption('minInnerIter', 1)
opt_test.optimizer.setOption('dynInnerIter', 1)
opt_test.optimizer.setOption('stopIters', 3)
opt_test.optimizer.setOption('w1', 0.99)
opt_test.optimizer.setOption('w2', 0.40)
opt_test.optimizer.setOption('c1', 2.0)
opt_test.optimizer.setOption('c2', 2.0)

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
opt_test.opt_prob.set_var('main_sole_z', value=0.3, lower=0.1, upper=0.4)

opt_test.opt_prob.set_covar('main_sole_z', 'main_sole_b')

# --------------------
# Run the optimization
opt_test.solve('mpirun -np 2 astra_r62_Linux_x86_64_OpenMPI_1.6.1')
# opt_test.solve('astra')
