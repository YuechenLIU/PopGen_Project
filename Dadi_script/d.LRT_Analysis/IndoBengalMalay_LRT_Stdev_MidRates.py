# !/usr/bin/python
# Filename: IndoBengalMalay_LRT_Stdev_MidRates.py
# Authorship: Yue-Chen Liu @ Peking University

"""A script to calculate the covariance of Indochina, Bengal, and Malay population."""

import dadi
import numpy

# define demographic models.

def IndoBengalMalay_TwoRates((nuIndo, nuBengal, nuMalay, T_IndoBengal, T_IndoMalay, mIndoBengalA, mBengalIndoA, mIndoBengalB, mBengalIndoB, mIndoMalay, mMalayIndo, mBengalMalay, mMalayBengal), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengalA, m21=mBengalIndoA)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)									# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoMalay, nu1=nuIndo, nu2=nuBengal, nu3=nuMalay, m12=mIndoBengalB, m21=mBengalIndoB, m13=mIndoMalay, m31=mMalayIndo, m23=mBengalMalay, m32=mMalayBengal)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))	# generate allele frequency spectrum.
	return fs
def IndoBengalMalay_MidRates((nuIndo, nuBengal, nuMalay, T_IndoBengal, T_IndoMalay, mIndoBengalA, mBengalIndoA, mIndoBengalB, mBengalIndoB, mIndoMalay, mMalayIndo), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengalA, m21=mBengalIndoA)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)									# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoMalay, nu1=nuIndo, nu2=nuBengal, nu3=nuMalay, m12=mIndoBengalB, m21=mBengalIndoB, m13=mIndoMalay, m31=mMalayIndo)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))									# generate allele frequency spectrum.
	return fs

# give the original data and bootstrap data.

data = dadi.Spectrum.from_file('../sfs_file/IndoBengalMalay.fs')
all_boot = [dadi.Spectrum.from_file('../sfs_file/IndoBengalMalay_bootstrap/{0:02d}.fs'.format(ii)) for ii in range(100)]
ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.

pts_l = [30,40,50]

# Choose a demographic model.

func_TwoRates = IndoBengalMalay_TwoRates
func_MidRates = IndoBengalMalay_MidRates

# Make the extrapolating version of our demographic model function.

func_ex_TwoRates = dadi.Numerics.make_extrap_log_func(func_TwoRates)
func_ex_MidRates = dadi.Numerics.make_extrap_log_func(func_MidRates)

# The followings are actual best-fit model parameters.

popt_TwoRates =	[1.91420276, 0.61752206, 0.10970641, 0.54291215, 0.04223735, 0.65461771, 3.50529324, 0.21500016, 0.00000001, 0.96906006, 6.28749694, 0.49288977, 0.73082084] # IndoBengalMalay_TwoRates.03.log
popt_MidRates =	[1.87310081, 0.63839652, 0.10247403, 0.54377007, 0.03651154, 0.90265292, 3.34938430, 0.00000140, 0.00000000, 0.00001117, 7.14118282] # IndoBengalMalay_MidRates.18.log

# Calculate the best-fit model AFS.

model_TwoRates = func_ex_TwoRates(popt_TwoRates, ns, pts_l)
model_MidRates = func_ex_MidRates(popt_MidRates, ns, pts_l)

# best-fit parameters found in multioptimizations.

ll_TwoRates = dadi.Inference.ll_multinom(model_TwoRates, data)
ll_MidRates = dadi.Inference.ll_multinom(model_MidRates, data)

# list of parameters for the complex model
# using the simple model best-fit params.

# MidRates:nuIndo, nuBengal, nuMalay, T_IndoBengal, T_IndoMalay, mIndoBengalA, mBengalIndoA, mIndoBengalB, mBengalIndoB, mIndoMalay, mMalayIndo
# TwoRates:nuIndo, nuBengal, nuMalay, T_IndoBengal, T_IndoMalay, mIndoBengalA, mBengalIndoA, mIndoBengalB, mBengalIndoB, mIndoMalay, mMalayIndo, mBengalMalay, mMalayBengal
p_lrt_TwoRates = [1.87310081, 0.63839652, 0.10247403, 0.54377007, 0.03651154, 0.90265292, 3.34938430, 0.00000140, 0.00000000, 0.00001117, 7.14118282, 0, 0]
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
