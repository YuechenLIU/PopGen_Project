# !/usr/bin/python
# Filename: IndoBengalAmur_noMig.py
# Authorship: Yue-Chen Liu @ Peking University

"""A python script to simulation the demographic model among Indochinese, Bengal, and Amur populations.
   This model only contain two population split events without migration among different groups.
   1st split between Bengal vs. Indochinese groups;
   2ed split between Amur vs. Indochinese groups.
"""

import dadi
# define a model.
def IndoBengalAmur_noMig((nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur), ns, pts):
        xx  = dadi.Numerics.default_grid(pts)           # specify the grid.
        phi = dadi.PhiManip.phi_1D(xx)                  # give the initial phi.
        phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)       # 1st split event.
        phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal)                # two population period.
        phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)                                       	# 2ed split event.
        phi = dadi.Integration.three_pops(phi, xx, T_IndoAmur, nu1=nuIndo, nu2=nuBengal, nu3=nuAmur)   	# three population period.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))        					# generate allele frequency spectrum.
        return fs

data = dadi.Spectrum.from_file('../sfs_file/IndoBengalAmur.fs')
ns = data.sample_sizes
pts_l = [30,40,50]
func = IndoBengalAmur_noMig

#parameters:  [nuIndo,	nuBengal,nuAmur,T_IndoBengal,T_IndoAmur]
upper_bound = [100,	100,	100,	10,	10]	# changed
lower_bound = [0.01,	0.01,	0.01,	0.00,	0.00]	# changed
p0          = [1.00,	1.00,	1.00,	0.10,	0.10]	# changed

func_ex = dadi.Numerics.make_extrap_log_func(func)
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, verbose=len(p0), maxiter=100)
print('Finshed optimization **************************************************')
print('Best-fit parameters: {0}'.format(popt))

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))
# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))
