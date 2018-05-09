# !/usr/bin/python
# Filename:IndoBengalMalyan_uncert.py
# Authorship: Yue-Chen Liu @ Peking University

"""A script to calculate covariance matrix of Indochina, Bengal, and Malayan population model."""

import dadi
import numpy

# clarify demographic models.

def IndoBengalMalay_noMig((nuIndo, nuBengal, nuMalay, T_IndoBengal, T_IndoMalay), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)								# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)									# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)							# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal)		# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)						# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoMalay, nu1=nuIndo, nu2=nuBengal, nu3=nuMalay)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))						# generate allele frequency spectrum.
	return fs
def IndoBengalMalay_13Mig((nuIndo, nuBengal, nuMalay, T_IndoBengal, T_IndoMalay, mIndoMalay, mMalayIndo), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)												# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)													# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)											# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal)						# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)										# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoMalay, nu1=nuIndo, nu2=nuBengal, nu3=nuMalay, m13=mIndoMalay, m31=mMalayIndo)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))										# generate allele frequency spectrum.
	return fs
def IndoBengalMalay_1213Mig((nuIndo, nuBengal, nuMalay, T_IndoBengal, T_IndoMalay, mIndoBengal, mBengalIndo, mIndoMalay, mMalayIndo), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)																	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)																		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)																# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengal, m21=mBengalIndo)							# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)															# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoMalay, nu1=nuIndo, nu2=nuBengal, nu3=nuMalay, m12=mIndoBengal, m21=mBengalIndo, m13=mIndoMalay, m31=mMalayIndo)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))															# generate allele frequency spectrum.
	return fs
def IndoBengalMalay_123Mig((nuIndo, nuBengal, nuMalay, T_IndoBengal, T_IndoMalay, mIndoBengal, mBengalIndo, mIndoMalay, mMalayIndo, mBengalMalay, mMalayBengal), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)		# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)			# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengal, m21=mBengalIndo)			# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)											# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoMalay, nu1=nuIndo, nu2=nuBengal, nu3=nuMalay, m12=mIndoBengal, m21=mBengalIndo, m13=mIndoMalay, m31=mMalayIndo, m23=mBengalMalay, m32=mMalayBengal)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))											# generate allele frequency spectrum.
	return fs
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
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)		# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuBengal, m12=mIndoBengalA, m21=mBengalIndoA)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)	# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoMalay, nu1=nuIndo, nu2=nuBengal, nu3=nuMalay, m12=mIndoBengalB, m21=mBengalIndoB, m13=mIndoMalay, m31=mMalayIndo)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))	# generate allele frequency spectrum.
	return fs

# give the original data and bootstrap data.

data = dadi.Spectrum.from_file('../sfs_file/IndoBengalMalay.fs')
all_boot = [dadi.Spectrum.from_file('../sfs_file/IndoBengalMalay_bootstrap/{0:02d}.fs'.format(ii)) for ii in range(100)]
ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.

pts_l = [30,40,50]

# Choose a demographic model.

func_noMig = IndoBengalMalay_noMig
func_13Mig = IndoBengalMalay_13Mig
func_1213Mig = IndoBengalMalay_1213Mig
func_123Mig = IndoBengalMalay_123Mig
func_TwoRates = IndoBengalMalay_TwoRates
func_MidRates = IndoBengalMalay_MidRates

# Make the extrapolating version of our demographic model function.

func_ex_noMig = dadi.Numerics.make_extrap_log_func(func_noMig)
func_ex_13Mig = dadi.Numerics.make_extrap_log_func(func_13Mig)
func_ex_1213Mig = dadi.Numerics.make_extrap_log_func(func_1213Mig)
func_ex_123Mig = dadi.Numerics.make_extrap_log_func(func_123Mig)
func_ex_TwoRates = dadi.Numerics.make_extrap_log_func(func_TwoRates)
func_ex_MidRates = dadi.Numerics.make_extrap_log_func(func_MidRates)

# The followings are actual best-fit model parameters

popt_noMig = [1.84822240, 0.76351764, 0.11148927, 0.14109691, 0.02765120] # IndoBengalMalyan_noMig.09.log
popt_13Mig = [1.86997707, 0.76421945, 0.09028367, 0.13613753, 0.03313133, 0.00000000, 8.44732062] # IndoBengalMalay_13Mig.07.log
popt_1213Mig = [1.92099292, 0.62247223, 0.09312841, 0.52925722, 0.03286610, 0.59696037, 2.44940901, 0.00002309, 8.06223989] # IndoBengalMalay_1213Mig.01.log
popt_123Mig = [1.92754173, 0.64907346, 0.08669454, 0.51976186, 0.03127931, 0.55378353, 2.16871987, 0.00010606, 7.53456650, 0.44024234, 1.05736785] # IndoBengalMalay_123Mig.18.log
popt_TwoRates = [1.91420276, 0.61752206, 0.10970641, 0.54291215, 0.04223735, 0.65461771, 3.50529324, 0.21500016, 0.00000001, 0.96906006, 6.28749694, 0.49288977, 0.73082084] # IndoBengalMalay_TwoRates.03.log
popt_MidRates = [1.87310081, 0.63839652, 0.10247403, 0.54377007, 0.03651154, 0.90265292, 3.34938430, 0.00000140, 0.00000000, 0.00001117, 7.14118282] # IndoBengalMalay_MidRates.18.log

# Calculate the best-fit model AFS.

model_noMig = func_ex_noMig(popt_noMig, ns, pts_l)
model_13Mig = func_ex_13Mig(popt_13Mig, ns, pts_l)
model_1213Mig = func_ex_1213Mig(popt_1213Mig, ns, pts_l)
model_123Mig = func_ex_123Mig(popt_123Mig, ns, pts_l)
model_TwoRates = func_ex_TwoRates(popt_TwoRates, ns, pts_l)
model_MidRates = func_ex_MidRates(popt_MidRates, ns, pts_l)

# Estimate parameter uncertainties using the Godambe Information Matrix.

GIM_noMig = dadi.Godambe.GIM_uncert(func_ex_noMig, pts_l, all_boot, popt_noMig, data, multinom=True, return_GIM=True)
GIM_13Mig = dadi.Godambe.GIM_uncert(func_ex_13Mig, pts_l, all_boot, popt_13Mig, data, multinom=True, return_GIM=True)
GIM_1213Mig = dadi.Godambe.GIM_uncert(func_ex_1213Mig, pts_l, all_boot, popt_1213Mig, data, multinom=True, return_GIM=True)
GIM_123Mig = dadi.Godambe.GIM_uncert(func_ex_123Mig, pts_l, all_boot, popt_123Mig, data, multinom=True, return_GIM=True)
GIM_TwoRates = dadi.Godambe.GIM_uncert(func_ex_TwoRates, pts_l, all_boot, popt_TwoRates, data, multinom=True, return_GIM=True)
GIM_MidRates = dadi.Godambe.GIM_uncert(func_ex_MidRates, pts_l, all_boot, popt_MidRates, data, multinom=True, return_GIM=True)

# Estimate the covarmatrix.

covarmatrix_noMig = numpy.linalg.inv(GIM_noMig[1])
covarmatrix_13Mig = numpy.linalg.inv(GIM_13Mig[1])
covarmatrix_1213Mig = numpy.linalg.inv(GIM_1213Mig[1])
covarmatrix_123Mig = numpy.linalg.inv(GIM_123Mig[1])
covarmatrix_TwoRates = numpy.linalg.inv(GIM_TwoRates[1])
covarmatrix_MidRates = numpy.linalg.inv(GIM_MidRates[1])

# Print the results.
print
print('Covariance matrix of noMig from GIM: {0}'.format(covarmatrix_noMig))
print
print('Covariance matrix of 13Mig from GIM: {0}'.format(covarmatrix_13Mig))
print
print('Covariance matrix of 1213Mig from GIM: {0}'.format(covarmatrix_1213Mig))
print
print('Covariance matrix of 123Mig from GIM: {0}'.format(covarmatrix_123Mig))
print
print('Covariance matrix of TwoRates from GIM: {0}'.format(covarmatrix_TwoRates))
print
print('Covariance matrix of MidRates from GIM: {0}'.format(covarmatrix_MidRates))
