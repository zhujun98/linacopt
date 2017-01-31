#!/usr/bin/python

"""
impact_test.py - A PYTHON script for beam matching optimization.

The optimization can bring down the mismatch factor from 1.41 to 1.17.
"""

import numpy as np
import sys
import os
sys.path.append(os.path.expanduser('~') + "/myscripts/linac_optimization/")

from linac_opt import LinacOpt

# -----------------------
# Instantiate the problem
matchin_opt = LinacOpt(path_name='impact_basic/',
                       input_file='ImpactT.in',
                       input_template='ImpactT.in.000',
                       particle_type='impact',
                       prob_name='matchin_opt')

# -----------------
# Set the optimizer
matchin_opt.set_optimizer('sdpen')
matchin_opt.optimizer.setOption('alfa_stop', 1e-3)
matchin_opt.optimizer.setOption('iprint', 0)
matchin_opt.optimizer.setOption('nf_max', 10000)

# --------------
# Add fit points
matchin_opt.fit_points.set_point('out', 'fort.107', q_norm=9.946e-17)

# -------------
# Add objective
def f1(fits):
    """Calculate the average mismatch factor."""
    mbetax = 0.04
    malphax = 1.6
    mbetay = 0.04
    malphay = 1.6
    mgammax = (1 + malphax**2)/mbetax
    mgammay = (1 + malphay**2)/mbetay

    betax = fits.out.betax
    alphax = fits.out.alphax
    betay = fits.out.betay
    alphay = fits.out.alphay
    gammax = (1 + alphax**2)/betax
    gammay = (1 + alphay**2)/betay

    bmag_x = 0.5*(betax*mgammax + gammax*mbetax - 2.0*alphax*malphax)
    bmag_y = 0.5*(betay*mgammay + gammay*mbetay - 2.0*alphay*malphay)

    mx = bmag_x + np.sqrt(bmag_x**2 - 1)
    my = bmag_y + np.sqrt(bmag_y**2 - 1)

    emitx = fits.out.emitx
    emity = fits.out.emity
    St = fits.out.St
    print("{:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f}".
          format((mx + my)/2, betax, alphax, betay, alphay, 1.0e6*emitx, 1.0e6*emity, 1.0e15*St))

    return (mx + my)/2

matchin_opt.opt_prob.set_obj('mismatch', func=f1)

# --------------
# Add constraint
def g1(fits):
    return fits.out.I_peak
matchin_opt.opt_prob.set_con('I_peak(A)', func=g1, lower=500)

def g2(fits):
    return fits.out.emitx*1e6
matchin_opt.opt_prob.set_con('emitx_um', func=g2, upper=1.0)

# ------------------------------------------------
# Add variables, co-variables and static-variables
matchin_opt.opt_prob.set_var('MQZM1_G', value= 0.5582, lower=-12.0, upper=12.0)
matchin_opt.opt_prob.set_var('MQZM2_G', value=-0.8456, lower=-12.0, upper=12.0)
matchin_opt.opt_prob.set_var('MQZM3_G', value= 0.4746, lower=-12.0, upper=12.0)
matchin_opt.opt_prob.set_var('MQZM4_G', value= 2.7909, lower=-12.0, upper=12.0)
matchin_opt.opt_prob.set_var('MQZM5_G', value=-5.2743, lower=-12.0, upper=12.0)
matchin_opt.opt_prob.set_var('MQZM6_G', value= 4.4487, lower=-12.0, upper=12.0)

matchin_opt.opt_prob.set_staticvar('n_col', 2)
matchin_opt.opt_prob.set_staticvar('n_row', 1)

# --------------------
# Run the optimization
# matchin_opt.solve('ImpactTv1.7linux')
matchin_opt.solve('mpirun -np 2 ImpactTv1.7linuxPara')
