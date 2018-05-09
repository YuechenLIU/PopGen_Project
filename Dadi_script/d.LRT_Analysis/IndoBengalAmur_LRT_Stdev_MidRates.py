# !/usr/bin/python
# Filename: IndoBengalAmur_LRT_Stdev_1213Mig.py
# Authorship: Yue-Chen Liu @ Peking University

"""A script to calculate the covariance of Indochina, Bengal, and Amur population."""

import dadi
import numpy

# define demographic models.

def IndoBengalAmur_TwoRates((nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengalA, mBengalIndoA, mIndoBengalB, mBengalIndoB, mIndoAmur, mAmurIndo, mBengalAmur, mAmurBengal), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengalA, m21=mBengalIndoA)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)									# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoAmur, nu1=nuIndo, nu2=nuBengal, nu3=nuAmur, m12=mIndoBengalB, m21=mBengalIndoB, m13=mIndoAmur, m31=mAmurIndo, m23=mBengalAmur, m32=mAmurBengal)
	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))	# generate allele frequency spectrum.
	return fs
def IndoBengalAmur_MidRates((nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengalA, mBengalIndoA, mIndoBengalB, mBengalIndoB, mIndoAmur, mAmurIndo), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengalA, m21=mBengalIndoA)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)									# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoAmur, nu1=nuIndo, nu2=nuBengal, nu3=nuAmur, m12=mIndoBengalB, m21=mBengalIndoB, m13=mIndoAmur, m31=mAmurIndo)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))									# generate allele frequency spectrum.
	return fs

# give the original data and bootstrap data.

data = dadi.Spectrum.from_file('../sfs_file/IndoBengalAmur.fs')
all_boot = [dadi.Spectrum.from_file('../sfs_file/IndoBengalAmur_bootstrap/{0:02d}.fs'.format(ii)) for ii in range(100)]
ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.

pts_l = [30,40,50]

# Choose a demographic model.

func_TwoRates = IndoBengalAmur_TwoRates
func_MidRates = IndoBengalAmur_MidRates

# Make the extrapolating version of our demographic model function.

func_ex_TwoRates = dadi.Numerics.make_extrap_log_func(func_TwoRates)
func_ex_MidRates = dadi.Numerics.make_extrap_log_func(func_MidRates)

# The followings are actual best-fit model parameters.

popt_TwoRates =	[1.92031158, 0.64592343, 0.09997775, 0.53120346, 0.07327762, 0.62862486, 4.37872560, 0.53004120, 0.00001072, 0.01855128, 0.09408632, 0.47844635, 0.29873964] # IndoBengalAmur_TwoRates.04.log
popt_MidRates = [1.89088687, 0.64940357, 0.09757292, 0.52694653, 0.06923142, 1.24695525, 4.98897616, 0.00078510, 0.00014790, 0.00006966, 0.18544790] # IndoBengalAmur_MidRates.03.log

# Calculate the best-fit model AFS.

model_TwoRates = func_ex_TwoRates(popt_TwoRates, ns, pts_l)
model_MidRates = func_ex_MidRates(popt_MidRates, ns, pts_l)

# best-fit parameters found in multioptimizations.

ll_TwoRates = dadi.Inference.ll_multinom(model_TwoRates, data)
ll_MidRates = dadi.Inference.ll_multinom(model_MidRates, data)

# list of parameters for the complex model
# using the simple model best-fit params.

# MidRates:nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengalA, mBengalIndoA, mIndoBengalB, mBengalIndoB, mIndoAmur, mAmurIndo
# TwoRates:nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengalA, mBengalIndoA, mIndoBengalB, mBengalIndoB, mIndoAmur, mAmurIndo, mBengalAmur, mAmurBengal
p_lrt_TwoRates = [1.89088687, 0.64940357, 0.09757292, 0.52694653, 0.06923142, 1.24695525, 4.98897616, 0.00078510, 0.00014790, 0.00006966, 0.18544790, 0, 0]
adjust_TwoRates = dadi.Godambe.LRT_adjust(func_ex_TwoRates, pts_l, all_boot, p_lrt_TwoRates, data, nested_indices=[11,12], multinom=True)
D_adj_TwoRates = adjust_TwoRates*2*(ll_TwoRates - ll_MidRates)
p_val_TwoRates = dadi.Godambe.sum_chi2_ppf(D_adj_TwoRates, weights=(0, 1))

# 1st output section

print
print('Maximum log composite likelihood of MidRates: {0}'.format(ll_MidRates))
print('Maximum log composite likelihood of TwoRates: {0}'.format(ll_TwoRates))

# 2ed output section

print
print('Adjusted D statistic of TwoRates and MidRates: {0:.4f}'.format(D_adj_TwoRates))
print('p-value for rejecting MidRates: {0:.4f}'.format(p_val_TwoRates))

# uncertain analysis section

uncerts_MidRates = dadi.Godambe.GIM_uncert(func_ex_MidRates, pts_l, all_boot, popt_MidRates, data, multinom=True)
uncerts_TwoRates = dadi.Godambe.GIM_uncert(func_ex_TwoRates, pts_l, all_boot, popt_TwoRates, data, multinom=True)

# 3rd output section

print
print('Estimated parameter standard deviations of MidRates from GIM: {0}'.format(uncerts_MidRates))
print('Estimated parameter standard deviations of TwoRates from GIM: {0}'.format(uncerts_TwoRates))
