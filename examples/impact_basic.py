#!/usr/bin/python

"""
This is a basic example showing how to optimize the beam size
in IMPACT-T with a local search optimizer.

It should end up with Sx = 0.0486 mm, Sy = 0.2998 mm.
"""
from linac_opt import LinacOpt


def f1(fits):
    """Define objective - transverse beam size"""
    print("Sx = {:8.4f} mm, Sy = {:8.4f} mm".
          format(1e3*fits.out.Sx, 1e3*fits.out.Sy))
    return 1e3*fits.out.Sx


def g1(fits):
    """Define constraint 1 - verticle beam size"""
    return 1e3*fits.out.Sy


# Instantiate the problem
matchin_opt = LinacOpt(path_name='impact_basic/',
                       input_file='ImpactT.in',
                       input_template='ImpactT.in.000',
                       particle_type='impact',
                       prob_name='matchin_opt')

# Set the optimizer
matchin_opt.set_optimizer('sdpen')
matchin_opt.optimizer.setOption('alfa_stop', 1e-2)
matchin_opt.optimizer.setOption('iprint', 0)
matchin_opt.optimizer.setOption('nf_max', 1000)

# Add fit points
matchin_opt.fit_points.set_point('out', 'fort.107', q_norm=9.946e-17)

# Add objective
matchin_opt.opt_prob.set_obj('mismatch', func=f1)

# Add constraints
matchin_opt.opt_prob.set_con('Sy(m)', func=g1, upper=0.3)
# Add variables
matchin_opt.opt_prob.set_var('MQZM1_G', value=0, lower=-12.0, upper=12.0)
matchin_opt.opt_prob.set_var('MQZM2_G', value=0, lower=-12.0, upper=12.0)

# Run the optimization
# matchin_opt.solve('ImpactTv1.7linux')
matchin_opt.opt_prob.set_staticvar('n_col', 2)
matchin_opt.opt_prob.set_staticvar('n_row', 1)
matchin_opt.solve('mpirun -np 2 ImpactTv1.7linuxPara')