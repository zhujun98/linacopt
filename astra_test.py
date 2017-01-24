#!/usr/bin/python
"""
astra_test.py - A PYTHON script for photoinjector optimization.

Space-charge is switched off for speed.
"""
import sys
import os
sys.path.append(os.path.expanduser('~') + "/myscripts/linac_optimization/")

from linac_opt import LinacOpt

# ----------------------------
# Instantiate the optimization

opt_test = LinacOpt(path_name='./astra_test',
                    input_file='injector.in',
                    input_template='injector.in.000',
                    particle_type='astra',
                    prob_name='opt_test',
                    restart=None,
                    run_once=False)

# ---------------------------------------------
# Set up the optimizer

opt_test.set_optimizer('alpso')

opt_test.optimizer.setOption('atol', 1e-2)
opt_test.optimizer.setOption('rtol', 1e-2)
opt_test.optimizer.setOption('etol', 1e-3)
opt_test.optimizer.setOption('itol', 1e-3)
opt_test.optimizer.setOption('SwarmSize', 10)  # For speed
opt_test.optimizer.setOption('maxOuterIter', 200)
opt_test.optimizer.setOption('maxInnerIter', 2)  # For speed
opt_test.optimizer.setOption('minInnerIter', 1)
opt_test.optimizer.setOption('dynInnerIter', 1)
opt_test.optimizer.setOption('stopIters', 2)
opt_test.optimizer.setOption('w1', 0.99)
opt_test.optimizer.setOption('w2', 0.40)
opt_test.optimizer.setOption('c1', 2.0)
opt_test.optimizer.setOption('c2', 2.0)
opt_test.optimizer.setOption('nf', 5)
opt_test.optimizer.setOption('ns', 15)

# --------------
# Add fit points

opt_test.fit_points.set_point('out', 'injector.0600.001', cut_tail=0.1)
opt_test.fit_points.set_point('slit', 'injector.0029.001', cut_halo=0.1)

# ------------
# Add sections

opt_test.sections.set_section('all')
opt_test.sections.set_section('tws2', z_lim=(3.7, 5.7))

# -------------
# Add objective

def f1(fits):
    print fits.out.St*1.e15
    return fits.out.St*1.e15
opt_test.opt_prob.set_obj('St_fs', f1)

# ---------------
# Add constraints

def g1(fits):
    return fits.out.n0
opt_test.opt_prob.set_con('n_pars', g1, equal=500, tol=0.0)

def g2(fits):
    return fits.slit.Sx*1.e3
opt_test.opt_prob.set_con('Sx_at_slit_mm', g2, upper=1.3)

def g3(fits, sections):
    return sections.all.Sx.max*1.e3
opt_test.opt_prob.set_con('max_Sx_mm', g3, upper=10.)

def g4(fits, sections):
    return sections.tws2.emitx.max*1.e6
opt_test.opt_prob.set_con('max_emitx_TWS2_um', g4, upper=0.3)

# ------------------------------------------------
# Add variables, co-variables and static-variables

opt_test.opt_prob.set_var('laser_spot', value=0.2, lower=0.1, upper=0.5)
opt_test.opt_prob.set_var('hc_sole_b', value=0.0, lower=0.0, upper=1.0)
opt_test.opt_prob.set_var('main_sole_b', value=0.0, lower=0.0, upper=0.4)
opt_test.opt_prob.set_var('tws1_sole_b', value=0.0, lower=0.0, upper=0.1)
opt_test.opt_prob.set_var('tws1_phase', value=0.0, lower=-90.0, upper=0.0)
opt_test.opt_prob.set_var('tws2_phase', value=0.0, lower=-90.0, upper=0.0)

opt_test.opt_prob.set_covar('tws2_sole_b', 'tws1_sole_b', slope=1.0, intercept=0.0)

opt_test.opt_prob.set_staticvar('hmax', value=0.02)
opt_test.opt_prob.set_staticvar('hmin', value=0.001)
opt_test.opt_prob.set_staticvar('nrad', value=16)
opt_test.opt_prob.set_staticvar('nlong', value=16)

# Code executive file name in your system
# run_code = 'ImpactTv1.7linux'
# run_code = 'mpirun -np 12 ImpactTv1.7linuxPara'
# run_code = 'astra'
run_code = 'mpirun -np 2 astra_r62_Linux_x86_64_OpenMPI_1.6.1'

# --------------------
# Run the optimization

opt_test.solve(run_code)

# *********************************************************************
# ************Refine using small step size and more grids**************
# *********************************************************************

# ----------------
# Change optimizer

opt_test.set_optimizer('sdpen')
opt_test.optimizer.setOption('alfa_stop', 1e-3)
opt_test.optimizer.setOption('iprint', 0)
opt_test.optimizer.setOption('nf_max', 200)

# ----------------
# Change objective

def f2(fits):
    print (fits.out.emitx + fits.out.emitx)*1.e6/2
    return (fits.out.emitx + fits.out.emitx)*1.e6/2
opt_test.opt_prob.set_obj('emit_um', f2)
opt_test.opt_prob.del_obj('St_fs')

# ---------------------------------------------------
# Change variables, co-variables and static-variables

opt_test.opt_prob.set_staticvar('hmax', value=0.01)
opt_test.opt_prob.set_staticvar('hmin', value=0.0005)
opt_test.opt_prob.set_staticvar('nrad', value=32)
opt_test.opt_prob.set_staticvar('nlong', value=32)

opt_test.opt_prob.set_staticvar('tws1_phase')
opt_test.opt_prob.set_staticvar('tws2_phase')

# ------------------
# Run the refinement

opt_test.solve(run_code)
