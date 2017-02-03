#!/usr/bin/python
"""
This example shows how to optimize the bunch length and emittance in
ASTRA using the velocity-bunching scheme. The global optimizer ALPSO
is used in the first step and the local search optimizer SDPEN is used
in the second one. Space-charge is switched off for speed.

Several advanced setups are also shown here.

The optimized result is 0.51 fs and 0.099 um. Make sure to set
SwarmSize == 80 and 'maxOuterIter' == 200.
"""
import sys
import os
sys.path.append(os.path.expanduser('~') + "/myscripts/linac_opt/")

from linac_opt import LinacOpt


# *********************************************************************
# ******************************Step 1*********************************
# *********************************************************************

# In the first step, the bunch length is going to be optimized with a
# global optimizer, i.e. the initial values of the variables take no
# effects.

# ----------------------------
# Instantiate the optimization
# - By setting restart to 1, the optimizer will try to read the saved
# result of the first step and directly start the second one.
# - If you are running a very unstable optimization, i.e. there are
# many crazy working points during the initial search, you may want
# to increase the value of max_fail (default is 20).
opt_test = LinacOpt(path_name='./astra_advanced',
                    input_file='injector.in',
                    input_template='injector.in.000',
                    particle_type='astra',
                    prob_name='opt_test',
                    restart=1,
                    max_fail=100,
                    run_once=False)

# ---------------------------------------------
# Set up the optimizer
opt_test.set_optimizer('alpso')

opt_test.optimizer.setOption('atol', 1e-2)
opt_test.optimizer.setOption('rtol', 1e-2)
opt_test.optimizer.setOption('SwarmSize', 10)
opt_test.optimizer.setOption('maxOuterIter', 3)
opt_test.optimizer.setOption('maxInnerIter', 3)  
opt_test.optimizer.setOption('minInnerIter', 1)
opt_test.optimizer.setOption('dynInnerIter', 1)
opt_test.optimizer.setOption('stopIters', 5)
opt_test.optimizer.setOption('w1', 0.99)
opt_test.optimizer.setOption('w2', 0.40)
opt_test.optimizer.setOption('c1', 2.0)
opt_test.optimizer.setOption('c2', 2.0)

# --------------
# Add fit points
opt_test.fit_points.set_point('out', 'injector.0600.001')
# The beam parameters at the 'slit' will be calculated after
# cutting 10% of the beam halo and then 10% of the beam tail.
opt_test.fit_points.set_point('slit', 'injector.0029.001',
                              cut_halo=0.1, cut_tail=0.1)

# ------------
# Add sections
opt_test.sections.set_section('all')
# 'tws2' is located between z = 3.7 ~ 5.7 m.
opt_test.sections.set_section('tws2', z_lim=(3.7, 5.7))


# -------------
# Add objective
def f1(fits):
    """Final bunch length"""
    print "rms bunch length: {:.2f} fs".format(fits.out.St*1.e15)
    return fits.out.St*1.e15
opt_test.opt_prob.set_obj('St_fs', f1)


# ---------------
# Add constraints
def g1(fits):
    """Final emittance"""
    return (fits.out.emitx + fits.out.emity)/2*1e6
opt_test.opt_prob.set_con('emitxy_um', g1, upper=0.30)


def g2(fits):
    """Beam size at the slit"""
    return fits.slit.Sx*1.e3
opt_test.opt_prob.set_con('Sx_at_slit_mm', g2, upper=2.0)


def g3(fits, sections):
    """Maximum beam size throughout the beamline"""
    return (sections.all.Sx.max + sections.tws2.Sy.max)/2*1.e3
opt_test.opt_prob.set_con('max_Sxy_mm', g3, upper=10.0)


def g4(fits, sections):
    """Average beta function inside TWS2"""
    return (sections.tws2.betax.ave + sections.tws2.betay.ave)/2
opt_test.opt_prob.set_con('ave_betaxy_TWS2_m', g4, equal=30.0, tol=20.0)

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
# The value of a co- variable can be dependant on several variables.
# opt_test.opt_prob.set_covar('tws2_sole_b', ('tws1_sole_b', 'tws1_sole_b'),
#                             slope=(1.0, 1.0), intercept=0.0)

opt_test.opt_prob.set_staticvar('charge', 0.010)  # 10 pC
opt_test.opt_prob.set_staticvar('laser_duration', 0.00018)  # 180 fs
opt_test.opt_prob.set_staticvar('hmax', 0.02)
opt_test.opt_prob.set_staticvar('hmin', 0.001)
opt_test.opt_prob.set_staticvar('nrad', 16)
opt_test.opt_prob.set_staticvar('nlong', 32)

# --------------------
# Run the optimization
# ASTRA will be killed after 10 seconds (default value is 1200 s).
# However, the optimization will continue. Both the shell command
# and the time_out will be inherited by the following optimizations.
# opt_test.solve('astra', time_out=10)
# The above input is equivalent to
opt_test.solve('timeout 10s astra injector.in >/dev/null', complete_shell=True)

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
