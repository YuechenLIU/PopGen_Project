# !/usr/bin/python
# Filename: IndoSumatranBengal_LRT_Stdev_noMig.py
# Authorship: Yue-Chen Liu @ Peking University

"""A script to calculate the covariance of Indochina, Bengal, and Amur population."""

import dadi
import numpy

# define demographic models.

def IndoSumatranBengal_noMig((nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)           # specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)                  # give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)       # 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoSumatran, nu1=nuIndo, nu2=nuSumatran)                    # two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)                                                       # 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuSumatran, nu3=nuBengal)      # three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))                                                        # generate allele frequency spectrum.
	return fs
def IndoSumatranBengal_13Mig((nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal, mIndoBengal, mBengalIndo), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)           # specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)                  # give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)       # 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoSumatran, nu1=nuIndo, nu2=nuSumatran)    # two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)                                       # 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuSumatran, nu3=nuBengal, m13=mIndoBengal, m31=mBengalIndo)   # three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))                                        # generate allele frequency spectrum.
	return fs
def IndoSumatranBengal_1213Mig((nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal, mIndoSumatran, mSumatranIndo, mIndoBengal, mBengalIndo), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)           # specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)                  # give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)       # 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoSumatran, nu1=nuIndo, nu2=nuSumatran, m12=mIndoSumatran, m21=mSumatranIndo)      # two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)                                                                               # 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuSumatran, nu3=nuBengal, m12=mIndoSumatran, m21=mSumatranIndo, m13=mIndoBengal, m31=mBengalIndo)      # three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))                                                                                        # generate allele frequency spectrum.
	return fs
def IndoSumatranBengal_123Mig((nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal, mIndoSumatran, mSumatranIndo, mIndoBengal, mBengalIndo, mSumatranBengal, mBengalSumatran), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)           # specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)                  # give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)       # 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoSumatran, nu1=nuIndo, nu2=nuSumatran, m12=mIndoSumatran, m21=mSumatranIndo)      # two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)                                                                               # 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuSumatran, nu3=nuBengal, m12=mIndoSumatran, m21=mSumatranIndo, m13=mIndoBengal, m31=mBengalIndo, m23=mSumatranBengal, m32=mBengalSumatran)
	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
	# generate allele frequency spectrum.
	return fs

# give the original data and bootstrap data.

data = dadi.Spectrum.from_file('../sfs_file/IndoSumatranBengal.fs')
all_boot = [dadi.Spectrum.from_file('../sfs_file/IndoSumatranBengal_bootstrap/{0:02d}.fs'.format(ii)) for ii in range(100)]
ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.

pts_l = [30,40,50]

# Choose a demographic model.

func_noMig = IndoSumatranBengal_noMig
func_13Mig = IndoSumatranBengal_13Mig
func_1213Mig = IndoSumatranBengal_1213Mig
func_123Mig = IndoSumatranBengal_123Mig

# Make the extrapolating version of our demographic model function.

func_ex_noMig = dadi.Numerics.make_extrap_log_func(func_noMig)
func_ex_13Mig = dadi.Numerics.make_extrap_log_func(func_13Mig)
func_ex_1213Mig = dadi.Numerics.make_extrap_log_func(func_1213Mig)
func_ex_123Mig = dadi.Numerics.make_extrap_log_func(func_123Mig)

# The followings are actual best-fit model parameters.

popt_noMig =	[1.97270365, 0.20137563, 0.71928058, 0.00075921, 0.16666017] # IndoSumatranBengal_noMig.07.log
popt_13Mig =	[1.81473904, 0.22290121, 0.72390819, 0.00000344, 0.18608391, 0.52416314, 0.20685316] # IndoSumatranBengal_13Mig.00.log
popt_1213Mig =	[1.69320015, 0.25786175, 0.70912252, 0.00000002, 0.27448429, 0.56412190, 0.31672346, 0.50463545, 1.43439920] # IndoSumatranBengal_1213Mig.01.log
popt_123Mig =	[1.64225874, 0.29668156, 0.67457777, 0.00001466, 0.40074425, 0.82522070, 0.32068534, 0.71792234, 1.41308053, 0.05999747, 0.58833900] # IndoSumatranBengal_123Mig.08.log

# Calculate the best-fit model AFS.

model_noMig = func_ex_noMig(popt_noMig, ns, pts_l)
model_13Mig = func_ex_13Mig(popt_13Mig, ns, pts_l)
model_1213Mig = func_ex_1213Mig(popt_1213Mig, ns, pts_l)
model_123Mig = func_ex_123Mig(popt_123Mig, ns, pts_l)

# best-fit parameters found in multioptimizations.

ll_noMig = dadi.Inference.ll_multinom(model_noMig, data)
ll_13Mig = dadi.Inference.ll_multinom(model_13Mig, data)
ll_1213Mig = dadi.Inference.ll_multinom(model_1213Mig, data)
ll_123Mig = dadi.Inference.ll_multinom(model_123Mig, data)

# list of parameters for the complex model
# using the simple model best-fit params.

# noMig: nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal
# 13Mig: nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal, mIndoBengal, mBengalIndo
p_lrt_13Mig = [1.97270365, 0.20137563, 0.71928058, 0.00075921, 0.16666017, 0, 0]
adjust_13Mig = dadi.Godambe.LRT_adjust(func_ex_13Mig, pts_l, all_boot, p_lrt_13Mig, data, nested_indices=[5,6], multinom=True)
D_adj_13Mig = adjust_13Mig*2*(ll_13Mig - ll_noMig)
p_val_13Mig = dadi.Godambe.sum_chi2_ppf(D_adj_13Mig, weights=(0, 1)) # one migration parameter plus 0.5 DOF.

# noMig:   nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal
# 1213Mig: nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal, mIndoSumatran, mSumatranIndo, mIndoBengal, mBengalIndo
p_lrt_1213Mig = [1.97270365, 0.20137563, 0.71928058, 0.00075921, 0.16666017, 0, 0, 0, 0]
adjust_1213Mig = dadi.Godambe.LRT_adjust(func_ex_1213Mig, pts_l, all_boot, p_lrt_1213Mig, data, nested_indices=[5,6,7,8], multinom=True)
D_adj_1213Mig = adjust_1213Mig*2*(ll_1213Mig - ll_noMig)
p_val_1213Mig = dadi.Godambe.sum_chi2_ppf(D_adj_1213Mig, weights=(0, 0, 1))

# noMig:   nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal
# 123Mig:  nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal, mIndoSumatran, mSumatranIndo, mIndoBengal, mBengalIndo, mSumatranBengal, mBengalSumatran
p_lrt_123Mig = [1.97270365, 0.20137563, 0.71928058, 0.00075921, 0.16666017, 0, 0, 0, 0, 0, 0]
adjust_123Mig = dadi.Godambe.LRT_adjust(func_ex_123Mig, pts_l, all_boot, p_lrt_123Mig, data, nested_indices=[5,6,7,8,9,10], multinom=True)
D_adj_123Mig = adjust_123Mig*2*(ll_123Mig - ll_noMig)
p_val_123Mig = dadi.Godambe.sum_chi2_ppf(D_adj_123Mig, weights=(0, 0, 0, 1))

# 1st output section

print
print('Maximum log composite likelihood of noMig: {0}'.format(ll_noMig))
print('Maximum log composite likelihood of 13Mig: {0}'.format(ll_13Mig))
print('Maximum log composite likelihood of 1213Mig: {0}'.format(ll_1213Mig))
print('Maximum log composite likelihood of 123Mig: {0}'.format(ll_123Mig))

# 2ed output section

print
print('Adjusted D statistic of 13Mig and noMig: {0:.4f}'.format(D_adj_13Mig))
print('p-value for rejecting noMig: {0:.4f}'.format(p_val_13Mig))
print
print('Adjusted D statistic of 1213Mig and noMig: {0:.4f}'.format(D_adj_1213Mig))
print('p-value for rejecting noMig: {0:.4f}'.format(p_val_1213Mig))
print
print('Adjusted D statistic of 123Mig and noMig: {0:.4f}'.format(D_adj_123Mig))
print('p-value for rejecting noMig: {0:.4f}'.format(p_val_123Mig))

# uncertain analysis section

uncerts_noMig = dadi.Godambe.GIM_uncert(func_ex_noMig, pts_l, all_boot, popt_noMig, data, multinom=True)
uncerts_13Mig = dadi.Godambe.GIM_uncert(func_ex_13Mig, pts_l, all_boot, popt_13Mig, data, multinom=True)
uncerts_1213Mig = dadi.Godambe.GIM_uncert(func_ex_1213Mig, pts_l, all_boot, popt_1213Mig, data, multinom=True)
uncerts_123Mig = dadi.Godambe.GIM_uncert(func_ex_123Mig, pts_l, all_boot, popt_123Mig, data, multinom=True)

# 3rd output section

print
print('Estimated parameter standard deviations of noMig from GIM: {0}'.format(uncerts_noMig))
print('Estimated parameter standard deviations of 13Mig from GIM: {0}'.format(uncerts_13Mig))
print('Estimated parameter standard deviations of 1213Mig from GIM: {0}'.format(uncerts_1213Mig))
print('Estimated parameter standard deviations of 123Mig from GIM: {0}'.format(uncerts_123Mig))
