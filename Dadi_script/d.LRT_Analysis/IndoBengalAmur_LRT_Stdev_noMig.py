# !/usr/bin/python
# Filename: IndoBengalAmur_LRT_Stdev_noMig.py
# Authorship: Yue-Chen Liu @ Peking University

"""A script to calculate the covariance of Indochina, Bengal, and Amur population."""

import dadi
import numpy

# define demographic models.

def IndoBengalAmur_noMig((nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)					# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoAmur, nu1=nuIndo, nu2=nuBengal, nu3=nuAmur)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))						# generate allele frequency spectrum.
	return fs
def IndoBengalAmur_12Mig((nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengal, mBengalIndo), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengal, m21=mBengalIndo)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)									# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoAmur, nu1=nuIndo, nu2=nuBengal, nu3=nuAmur, m12=mIndoBengal, m21=mBengalIndo)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))										# generate allele frequency spectrum.
	return fs
def IndoBengalAmur_1213Mig((nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengal, mBengalIndo, mIndoAmur, mAmurIndo), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengal, m21=mBengalIndo)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)									# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoAmur, nu1=nuIndo, nu2=nuBengal, nu3=nuAmur, m12=mIndoBengal, m21=mBengalIndo, m13=mIndoAmur, m31=mAmurIndo)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))														# generate allele frequency spectrum.
	return fs
def IndoBengalAmur_123Mig((nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengal, mBengalIndo, mIndoAmur, mAmurIndo, mBengalAmur, mAmurBengal), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)		# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)			# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengal, m21=mBengalIndo)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)									# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoAmur, nu1=nuIndo, nu2=nuBengal, nu3=nuAmur, m12=mIndoBengal, m21=mBengalIndo, m13=mIndoAmur, m31=mAmurIndo, m23=mBengalAmur, m32=mAmurBengal)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))	# generate allele frequency spectrum.
	return fs
def IndoBengalAmur_TwoRates((nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengalA, mBengalIndoA, mIndoBengalB, mBengalIndoB, mIndoAmur, mAmurIndo, mBengalAmur, mAmurBengal), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengalA, m21=mBengalIndoA)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)									# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoAmur, nu1=nuIndo, nu2=nuBengal, nu3=nuAmur, m12=mIndoBengalB, m21=mBengalIndoB, m13=mIndoAmur, m31=mAmurIndo, m23=mBengalAmur, m32=mAmurBengal)	# three population period.
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

func_noMig = IndoBengalAmur_noMig
func_12Mig = IndoBengalAmur_12Mig
func_1213Mig = IndoBengalAmur_1213Mig
func_123Mig = IndoBengalAmur_123Mig
func_TwoRates = IndoBengalAmur_TwoRates
func_MidRates = IndoBengalAmur_MidRates

# Make the extrapolating version of our demographic model function.

func_ex_noMig = dadi.Numerics.make_extrap_log_func(func_noMig)
func_ex_12Mig = dadi.Numerics.make_extrap_log_func(func_12Mig)
func_ex_1213Mig = dadi.Numerics.make_extrap_log_func(func_1213Mig)
func_ex_123Mig = dadi.Numerics.make_extrap_log_func(func_123Mig)
func_ex_TwoRates = dadi.Numerics.make_extrap_log_func(func_TwoRates)
func_ex_MidRates = dadi.Numerics.make_extrap_log_func(func_MidRates)

# The followings are actual best-fit model parameters.

popt_noMig =	[1.91444067, 0.79876312, 0.09318636, 0.11430078, 0.06548706] # IndoBengalAmur_noMig.17.log
popt_12Mig =	[1.95667364, 0.59896642, 0.09989642, 0.50330380, 0.06868595, 0.50227995, 2.60881723] # IndoBengalAmur_12Mig.00.log
popt_1213Mig =	[1.96633697, 0.59204999, 0.09947258, 0.50197729, 0.07074859, 0.48509181, 2.66006880, 0.00034916, 0.27342419] # IndoBengalAmur_1213Mig.11.log
popt_123Mig =	[1.93405202, 0.65427995, 0.09620943, 0.48653942, 0.07045946, 0.55807477, 2.00647798, 0.00000015, 0.05170904, 0.58705115, 0.37247137] # IndoBengalAmur_123Mig.16.log
popt_TwoRates =	[1.92031158, 0.64592343, 0.09997775, 0.53120346, 0.07327762, 0.62862486, 4.37872560, 0.53004120, 0.00001072, 0.01855128, 0.09408632, 0.47844635, 0.29873964] # IndoBengalAmur_TwoRates.04.log
popt_MidRates = [1.89088687, 0.64940357, 0.09757292, 0.52694653, 0.06923142, 1.24695525, 4.98897616, 0.00078510, 0.00014790, 0.00006966, 0.18544790] # IndoBengalAmur_MidRates.03.log

# Calculate the best-fit model AFS.

model_noMig = func_ex_noMig(popt_noMig, ns, pts_l)
model_12Mig = func_ex_12Mig(popt_12Mig, ns, pts_l)
model_1213Mig = func_ex_1213Mig(popt_1213Mig, ns, pts_l)
model_123Mig = func_ex_123Mig(popt_123Mig, ns, pts_l)
model_TwoRates = func_ex_TwoRates(popt_TwoRates, ns, pts_l)
model_MidRates = func_ex_MidRates(popt_MidRates, ns, pts_l)

# best-fit parameters found in multioptimizations.

ll_noMig = dadi.Inference.ll_multinom(model_noMig, data)
ll_12Mig = dadi.Inference.ll_multinom(model_12Mig, data)
ll_1213Mig = dadi.Inference.ll_multinom(model_1213Mig, data)
ll_123Mig = dadi.Inference.ll_multinom(model_123Mig, data)
ll_TwoRates = dadi.Inference.ll_multinom(model_TwoRates, data)
ll_MidRates = dadi.Inference.ll_multinom(model_MidRates, data)

# list of parameters for the complex model
# using the simple model best-fit params.

# noMig: nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur
# 12Mig: nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengal, mBengalIndo
p_lrt_12Mig = [1.91444067, 0.79876312, 0.09318636, 0.11430078, 0.06548706, 0, 0]
adjust_12Mig = dadi.Godambe.LRT_adjust(func_ex_12Mig, pts_l, all_boot, p_lrt_12Mig, data, nested_indices=[5,6], multinom=True)
D_adj_12Mig = adjust_12Mig*2*(ll_12Mig - ll_noMig)
p_val_12Mig = dadi.Godambe.sum_chi2_ppf(D_adj_12Mig, weights=(0, 1)) # one migration parameter plus 0.5 DOF.

# noMig:   nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur
# 1213Mig: nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengal, mBengalIndo, mIndoAmur, mAmurIndo
p_lrt_1213Mig = [1.91444067, 0.79876312, 0.09318636, 0.11430078, 0.06548706, 0, 0, 0, 0]
adjust_1213Mig = dadi.Godambe.LRT_adjust(func_ex_1213Mig, pts_l, all_boot, p_lrt_1213Mig, data, nested_indices=[5,6,7,8], multinom=True)
D_adj_1213Mig = adjust_1213Mig*2*(ll_1213Mig - ll_noMig)
p_val_1213Mig = dadi.Godambe.sum_chi2_ppf(D_adj_1213Mig, weights=(0, 0, 1))

# noMig:   nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur
# 123Mig:  nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengal, mBengalIndo, mIndoAmur, mAmurIndo, mBengalAmur, mAmurBengal
p_lrt_123Mig = [1.91444067, 0.79876312, 0.09318636, 0.11430078, 0.06548706, 0, 0, 0, 0, 0, 0]
adjust_123Mig = dadi.Godambe.LRT_adjust(func_ex_123Mig, pts_l, all_boot, p_lrt_123Mig, data, nested_indices=[5,6,7,8,9,10], multinom=True)
D_adj_123Mig = adjust_123Mig*2*(ll_123Mig - ll_noMig)
p_val_123Mig = dadi.Godambe.sum_chi2_ppf(D_adj_123Mig, weights=(0, 0, 0, 1))

# noMig:   nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur
# MidRates:nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengalA, mBengalIndoA, mIndoBengalB, mBengalIndoB, mIndoAmur, mAmurIndo
p_lrt_MidRates = [1.91444067, 0.79876312, 0.09318636, 0.11430078, 0.06548706, 0, 0, 0, 0, 0, 0]
adjust_MidRates = dadi.Godambe.LRT_adjust(func_ex_MidRates, pts_l, all_boot, p_lrt_MidRates, data, nested_indices=[5,6,7,8,9,10], multinom=True)
D_adj_MidRates = adjust_MidRates*2*(ll_MidRates - ll_noMig)
p_val_MidRates = dadi.Godambe.sum_chi2_ppf(D_adj_MidRates, weights=(0, 0, 0, 1))

# noMig:   nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur
# TwoRates:nuIndo, nuBengal, nuAmur, T_IndoBengal, T_IndoAmur, mIndoBengalA, mBengalIndoB, mIndoBengalA, mBengalIndoB, mIndoAmur, mAmurIndo, mBengalAmur, mAmurBengal
p_lrt_TwoRates = [1.91444067, 0.79876312, 0.09318636, 0.11430078, 0.06548706, 0, 0, 0, 0, 0, 0, 0, 0]
adjust_TwoRates = dadi.Godambe.LRT_adjust(func_ex_TwoRates, pts_l, all_boot, p_lrt_TwoRates, data, nested_indices=[5,6,7,8,9,10,11,12], multinom=True)
D_adj_TwoRates = adjust_TwoRates*2*(ll_TwoRates - ll_noMig)
p_val_TwoRates = dadi.Godambe.sum_chi2_ppf(D_adj_TwoRates, weights=(0, 0, 0, 0, 1))

# 1st output section

print
print('Maximum log composite likelihood of noMig: {0}'.format(ll_noMig))
print('Maximum log composite likelihood of 12Mig: {0}'.format(ll_12Mig))
print('Maximum log composite likelihood of 1213Mig: {0}'.format(ll_1213Mig))
print('Maximum log composite likelihood of 123Mig: {0}'.format(ll_123Mig))
print('Maximum log composite likelihood of MidRates: {0}'.format(ll_MidRates))
print('Maximum log composite likelihood of TwoRates: {0}'.format(ll_TwoRates))

# 2ed output section

print
print('Adjusted D statistic of 12Mig and noMig: {0:.4f}'.format(D_adj_12Mig))
print('p-value for rejecting noMig: {0:.4f}'.format(p_val_12Mig))
print
print('Adjusted D statistic of 1213Mig and noMig: {0:.4f}'.format(D_adj_1213Mig))
print('p-value for rejecting noMig: {0:.4f}'.format(p_val_1213Mig))
print
print('Adjusted D statistic of 123Mig and noMig: {0:.4f}'.format(D_adj_123Mig))
print('p-value for rejecting noMig: {0:.4f}'.format(p_val_123Mig))
print
print('Adjusted D statistic of MidRates and noMig: {0:.4f}'.format(D_adj_MidRates))
print('p-value for rejecting noMig: {0:.4f}'.format(p_val_MidRates))
print
print('Adjusted D statistic of TwoRates and noMig: {0:.4f}'.format(D_adj_TwoRates))
print('p-value for rejecting noMig: {0:.4f}'.format(p_val_TwoRates))

# uncertain analysis section

uncerts_noMig = dadi.Godambe.GIM_uncert(func_ex_noMig, pts_l, all_boot, popt_noMig, data, multinom=True)
uncerts_12Mig = dadi.Godambe.GIM_uncert(func_ex_12Mig, pts_l, all_boot, popt_12Mig, data, multinom=True)
uncerts_1213Mig = dadi.Godambe.GIM_uncert(func_ex_1213Mig, pts_l, all_boot, popt_1213Mig, data, multinom=True)
uncerts_123Mig = dadi.Godambe.GIM_uncert(func_ex_123Mig, pts_l, all_boot, popt_123Mig, data, multinom=True)
uncerts_MidRates = dadi.Godambe.GIM_uncert(func_ex_MidRates, pts_l, all_boot, popt_MidRates, data, multinom=True)
uncerts_TwoRates = dadi.Godambe.GIM_uncert(func_ex_TwoRates, pts_l, all_boot, popt_TwoRates, data, multinom=True)

# 3rd output section

print
print('Estimated parameter standard deviations of noMig from GIM: {0}'.format(uncerts_noMig))
print('Estimated parameter standard deviations of 12Mig from GIM: {0}'.format(uncerts_12Mig))
print('Estimated parameter standard deviations of 1213Mig from GIM: {0}'.format(uncerts_1213Mig))
print('Estimated parameter standard deviations of 123Mig from GIM: {0}'.format(uncerts_123Mig))
print('Estimated parameter standard deviations of MidRates from GIM: {0}'.format(uncerts_MidRates))
print('Estimated parameter standard deviations of TwoRates from GIM: {0}'.format(uncerts_TwoRates))
