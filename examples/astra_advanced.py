#!/usr/bin/python
"""
astra_test.py - A PYTHON script for photoinjector optimization.

Space-charge is switched off for speed.
"""
import sys
import os
sys.path.append(os.path.expanduser('~') + "/myscripts/linac_optimization/")

from linac_opt import LinacOpt


# *********************************************************************
# ******************************Step 1*********************************
# *********************************************************************

# In the first step, the bunch length is going to be optimized with a
# global optimizer, i.e. the initial values of the variables have no
# effects.

# ----------------------------
# Instantiate the optimization
opt_test = LinacOpt(path_name='./astra_advanced',
                    input_file='injector.in',
                    input_template='injector.in.000',
                    particle_type='astra',
                    prob_name='opt_test',
                    restart=None)

# ---------------------------------------------
# Set up the optimizer
opt_test.set_optimizer('alpso')

opt_test.optimizer.setOption('atol', 1e-1)
opt_test.optimizer.setOption('rtol', 1e-2)
opt_test.optimizer.setOption('SwarmSize', 40)  
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
opt_test.fit_points.set_point('out', 'injector.0600.001', cut_tail=0.1)
opt_test.fit_points.set_point('slit', 'injector.0029.001', cut_halo=0.1)

# ------------
# Add sections
# opt_test.sections.set_section('all')
# opt_test.sections.set_section('tws2', z_lim=(3.7, 5.7))

# -------------
# Add objective
def f1(fits):
    print "rms bunch length: {:.2f} fs".format(fits.out.St*1.e15)
    return fits.out.St*1.e15
opt_test.opt_prob.set_obj('St_fs', f1)

# ---------------
# Add constraints
# def g1(fits):
#     return (fits.out.emitx + fits.out.emity)/2*1e6
# opt_test.opt_prob.set_con('emitxy_um', g1, upper=0.20)
#
# def g2(fits):
#     return fits.slit.Sx*1.e3
# opt_test.opt_prob.set_con('Sx_at_slit_mm', g2, upper=2.0)
#
# def g3(fits, sections):
#     return (sections.all.Sx.max + sections.tws2.Sy.max)/2*1.e3
# opt_test.opt_prob.set_con('max_Sxy_mm', g3, upper=10.0)
#
# def g4(fits, sections):
#     return (sections.tws2.betax.ave + sections.tws2.betay.ave)/2
# opt_test.opt_prob.set_con('ave_betaxy_TWS2_m', g4, equal=20.0, tol=10.0)

# ------------------------------------------------
# Add variables, co-variables and static-variables
opt_test.opt_prob.set_var('gun_phase', value=0.0, lower=-30.0, upper=30.0)
opt_test.opt_prob.set_var('laser_spot', value=0.2, lower=0.1, upper=0.5)
opt_test.opt_prob.set_var('hc_sole_b', value=0.0, lower=0.0, upper=0.6)
opt_test.opt_prob.set_var('main_sole_b', value=0.0, lower=0.0, upper=0.4)
opt_test.opt_prob.set_var('tws1_sole_b', value=0.0, lower=0.0, upper=0.1)
opt_test.opt_prob.set_var('tws1_phase', value=0.0, lower=-90.0, upper=0.0)
opt_test.opt_prob.set_var('tws2_phase', value=0.0, lower=-90.0, upper=0.0)

opt_test.opt_prob.set_covar('tws2_sole_b', 'tws1_sole_b', slope=1.0, intercept=0.0)

opt_test.opt_prob.set_staticvar('charge', 0.010)  # 10 pC
opt_test.opt_prob.set_staticvar('laser_duration', 0.00018)  # 180 fs
opt_test.opt_prob.set_staticvar('hmax', 0.02)
opt_test.opt_prob.set_staticvar('hmin', 0.001)
opt_test.opt_prob.set_staticvar('nrad', 16)
opt_test.opt_prob.set_staticvar('nlong', 32)

# --------------------
# Run the optimization
# ASTRA will be killed after 10 seconds. However, the optimization
# will continue. Both the shell command and the max_duration will
# be inherited by the following optimizations.
opt_test.solve('astra', max_duration=10)
# opt_test.solve('mpirun -np 2 astra_r62_Linux_x86_64_OpenMPI_1.6.1')

# *********************************************************************
# ******************************Step 2*********************************
# *********************************************************************

# In the second step, the optimizer is set to the local search optimizer
# SDPEN and the objective is changed to emittance. The phases of the gun
# and the TWSs are also removed from the variable list since we do not
# want to change the bunch length.

# ----------------
# Change optimizer
opt_test.set_optimizer('sdpen')
opt_test.optimizer.setOption('alfa_stop', 1e-2)
opt_test.optimizer.setOption('iprint', 0)
opt_test.optimizer.setOption('nf_max', 5000)

# ----------------
# Change objective
def f2(fits):
    print (fits.out.emitx + fits.out.emitx)*1.e6/2
    return (fits.out.emitx + fits.out.emitx)*1.e6/2
opt_test.opt_prob.set_obj('emit_um', f2)
opt_test.opt_prob.del_obj('St_fs')

# ---------------------------------------------------
# Change variables, co-variables and static-variables
opt_test.opt_prob.set_staticvar('hmax', 0.01)
opt_test.opt_prob.set_staticvar('hmin', 0.0005)
opt_test.opt_prob.set_staticvar('nrad', 32)
opt_test.opt_prob.set_staticvar('nlong', 64)

opt_test.opt_prob.set_staticvar('gun_phase')
opt_test.opt_prob.set_staticvar('tws1_phase')
opt_test.opt_prob.set_staticvar('tws2_phase')

# ------------------
# Run the refinement
opt_test.solve()