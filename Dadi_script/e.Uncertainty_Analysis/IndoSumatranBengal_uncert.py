# !/usr/bin/python
# Filename: IndoSumatranBengal_uncert.py 
# Authorship: Yue-Chen Liu @ Peking University

"""A script to calculate the covariance of Indochina, Sumatran and Bengal population models."""
import dadi
import numpy

# define demographic models.

def IndoSumatranBengal_noMig((nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoSumatran, nu1=nuIndo, nu2=nuSumatran)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)					# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuSumatran, nu3=nuBengal)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))							# generate allele frequency spectrum.
	return fs
def IndoSumatranBengal_13Mig((nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal, mIndoBengal, mBengalIndo), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		#  give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoSumatran, nu1=nuIndo, nu2=nuSumatran)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)					# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuSumatran, nu3=nuBengal, m13=mIndoBengal, m31=mBengalIndo)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))					# generate allele frequency spectrum.
	return fs
def IndoSumatranBengal_1213Mig((nuIndo, nuSumatran, nuBengal, T_IndoSumatran, T_IndoBengal, mIndoSumatran, mSumatranIndo, mIndoBengal, mBengalIndo), ns, pts):
	xx  = dadi.Numerics.default_grid(pts)	# specify the grid.
	phi = dadi.PhiManip.phi_1D(xx)		# give the initial phi.
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)	# 1st split event.
	phi = dadi.Integration.two_pops(phi, xx, T_IndoSumatran, nu1=nuIndo, nu2=nuSumatran, m12=mIndoSumatran, m21=mSumatranIndo)	# two population period.
	phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)										# 2ed split event.
	phi = dadi.Integration.three_pops(phi, xx, T_IndoBengal, nu1=nuIndo, nu2=nuSumatran, nu3=nuBengal, m12=mIndoSumatran, m21=mSumatranIndo, m13=mIndoBengal, m31=mBengalIndo)	# three population period.
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))	# generate allele frequency spectrum.
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

# Estimate parameter uncertainties using the Godambe Information Matrix.

GIM_noMig = dadi.Godambe.GIM_uncert(func_ex_noMig, pts_l, all_boot, popt_noMig, data, multinom=True, return_GIM=True)
GIM_13Mig = dadi.Godambe.GIM_uncert(func_ex_13Mig, pts_l, all_boot, popt_13Mig, data, multinom=True, return_GIM=True)
GIM_1213Mig = dadi.Godambe.GIM_uncert(func_ex_1213Mig, pts_l, all_boot, popt_1213Mig, data, multinom=True, return_GIM=True)
GIM_123Mig = dadi.Godambe.GIM_uncert(func_ex_123Mig, pts_l, all_boot, popt_123Mig, data, multinom=True, return_GIM=True)

# Estimate the covarmatrix.

covarmatrix_noMig = numpy.linalg.inv(GIM_noMig[1])
covarmatrix_13Mig = numpy.linalg.inv(GIM_13Mig[1])
covarmatrix_1213Mig = numpy.linalg.inv(GIM_1213Mig[1])
covarmatrix_123Mig = numpy.linalg.inv(GIM_123Mig[1])

# Print the results.
print
print('Covariance matrix of noMig from GIM: {0}'.format(covarmatrix_noMig))
print
print('Covariance matrix of 13Mig from GIM: {0}'.format(covarmatrix_13Mig))
print
print('Covariance matrix of 1213Mig from GIM: {0}'.format(covarmatrix_1213Mig))
print
print('Covariance matrix of 123Mig from GIM: {0}'.format(covarmatrix_123Mig))
