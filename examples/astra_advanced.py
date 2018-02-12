#!/usr/bin/python
"""
This example shows how to optimize the bunch length and emittance in
ASTRA using the velocity-bunching scheme. The global optimizer ALPSO
is used in the first step and the local search optimizer SDPEN is used
in the second one. Space-charge is switched off for speed.

Several advanced setups are also shown here.
"""
from linacopt import LinacOpt


# ******************************Step 1*********************************
# In the first step, the bunch length is going to be optimized with a
# global optimizer, i.e. the initial values of the variables take no
# effects.
# *********************************************************************

def f1(fits):
    """Define objective - final bunch length"""
    print("rms bunch length: {:.2f} fs".format(fits.out.St*1.e15))
    return fits.out.St*1.e15


def g1(fits):
    """Define constraint 1 - final emittance"""
    return (fits.out.emitx + fits.out.emity)/2*1e6


def g2(fits):
    """Define constraint 2 - beam size at the slit"""
    return fits.slit.Sx*1.e3


def g3(fits, sections):
    """Define constraint 3 - maximum beam size throughout the beamline

    The second argument 'sections' refers to a section in the beamline,
    and it can be the whole beam line. The statistic parameters in one
    'section' can be used for optimization. See beam_evolutions.py.

    Note: fits (fit points) must be the first argument!
    """
    return (sections.all.Sx.max + sections.all.Sy.max)/2*1.e3


def g4(fits, sections):
    """Define constraint 4 - average beta function inside TWS2"""
    return (sections.tws2.betax.ave + sections.tws2.betay.ave)/2


# Instantiate the optimization
opt_test = LinacOpt(path_name='./astra_advanced',
                    input_file='injector.in',
                    input_template='injector.in.000',
                    particle_type='astra',
                    restart=None,
                    max_fail=100)

# Set up the optimizer
opt_test.set_optimizer('alpso', pll_type='SPM')

opt_test.optimizer.setOption('atol', 1e-2)
opt_test.optimizer.setOption('rtol', 1e-2)
opt_test.optimizer.setOption('SwarmSize', 80)
opt_test.optimizer.setOption('maxOuterIter', 30)
opt_test.optimizer.setOption('maxInnerIter', 3)  
opt_test.optimizer.setOption('minInnerIter', 1)
opt_test.optimizer.setOption('dynInnerIter', 1)
opt_test.optimizer.setOption('stopIters', 5)
opt_test.optimizer.setOption('w1', 0.99)
opt_test.optimizer.setOption('w2', 0.40)
opt_test.optimizer.setOption('c1', 2.0)
opt_test.optimizer.setOption('c2', 2.0)

# fit points
opt_test.fit_points.set_point('out', 'injector.0620.001')
opt_test.fit_points.set_point('slit', 'injector.0029.001',
                              cut_halo=0.1, cut_tail=0.1, rotate=20.0)

# objective
opt_test.opt_prob.set_obj('St_fs', f1)

# sections
opt_test.sections.set_section('all')
opt_test.sections.set_section('tws2', z_lim=(4.0, 6.0))

# constraints
opt_test.opt_prob.set_con('emitxy_um', g1, upper=0.30)
opt_test.opt_prob.set_con('Sx_at_slit_mm', g2, upper=2.0)
opt_test.opt_prob.set_con('max_Sxy_mm', g3, upper=10.0)
opt_test.opt_prob.set_con('ave_betaxy_TWS2_m', g4, equal=30.0, tol=20.0)

# variables (global optimization, no need to set initial values)
opt_test.opt_prob.set_var('laser_spot', lower=0.1, upper=0.5)
opt_test.opt_prob.set_var('laser_duration', lower=0.0002, upper=0.001)
opt_test.opt_prob.set_var('main_sole_b', lower=0.0, upper=0.4)
opt_test.opt_prob.set_var('tws1_sole_b', lower=0.0, upper=0.1)
opt_test.opt_prob.set_var('gun_phase', lower=-30.0, upper=30.0)
opt_test.opt_prob.set_var('tws1_phase', lower=-90.0, upper=0.0)
opt_test.opt_prob.set_var('tws2_phase', lower=-90.0, upper=0.0)

# co-variables
opt_test.opt_prob.set_covar('tws2_sole_b', 'tws1_sole_b', 1.0, 0.0)

# static-variables
opt_test.opt_prob.set_staticvar('gun_gradient', 110)
opt_test.opt_prob.set_staticvar('tws1_gradient', 30.0)
opt_test.opt_prob.set_staticvar('tws2_gradient', 30.0)
opt_test.opt_prob.set_staticvar('charge', 0.010)
opt_test.opt_prob.set_staticvar('hmax', 0.005)
opt_test.opt_prob.set_staticvar('hmin', 0.00005)
opt_test.opt_prob.set_staticvar('nrad', 16)
opt_test.opt_prob.set_staticvar('nlong', 32)

# Run the optimization
#
# The parallel version of ASTRA will get stuck sometimes when a lot
# of particles are lost. Hence, ASTRA will be killed after 5 seconds
# (assuming usually one iteration takes 3 seconds). However, the
# optimization will continue. Both the shell command and the time_out
# will be inherited by the following optimizations.
# opt_test.solve('astra')
opt_test.solve('mpirun -np 2 astra_r62_Linux_x86_64_OpenMPI_1.6.1')


# ******************************Step 2*********************************


def f2(fits):
    """Define new objective"""
    print("Transverse emittance: {}".format((fits.out.emitx + fits.out.emitx)*1.0e6/2))
    return (fits.out.emitx + fits.out.emitx)*1.e6/2

# Change optimizer
#
# In the second step, the optimizer is set to the local search optimizer
# SDPEN and the objective is changed to emittance. The phases and
# amplitudes of the gun and the TWSs are also removed from the variable
# list since we do not want to change the bunch length a lot.
opt_test.set_optimizer('sdpen')
opt_test.optimizer.setOption('alfa_stop', 1e-3)

# Change objective
opt_test.opt_prob.set_obj('emit_um', f2)
opt_test.opt_prob.del_obj('St_fs')

# Add another constraint
opt_test.opt_prob.set_con('St_fs', f1, upper=0.1)

# Change variables to static variables
# Remove variables to which the bunch length is sensitive. (You can also skip this step)
opt_test.opt_prob.set_staticvar('gun_phase')
opt_test.opt_prob.set_staticvar('tws1_phase')
opt_test.opt_prob.set_staticvar('tws2_phase')

# Run the local search optimization
opt_test.solve()
