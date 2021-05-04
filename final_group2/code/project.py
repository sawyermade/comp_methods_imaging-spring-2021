import numpy as np, scipy.io, sys, os
import matplotlib.pyplot as plt

def define_psf(U, V, slope):
	# Vars
	x = np.linspace(-1, 1, 2 * U).astype(float)
	y = np.linspace(-1, 1, 2 * U).astype(float)
	z = np.linspace(0, 2, 2 * V).astype(float)
	# gz, gy, gx = [np.swapaxes(i, 0, 1) for i in np.meshgrid(z, y, x)]
	gz, gy, gx = np.meshgrid(z, y, x, indexing='ij')
	print(f'x, y, z shape: {x.shape}, {y.shape}, {z.shape}')
	print(f'gx, gy, gz shape: {gx.shape}, {gy.shape}, {gz.shape}')
	print(f'gx type: {type(gx)}')

	# Define PSF
	psf = np.abs((4 * slope)**2 * (gx**2 + gy**2) - gz).astype(float)
	print(f'psf shape: {psf.shape}')
	psf = np.tile(psf.min(0), [2 * V, 1, 1])
	print(f'psf tiled shape: {psf.shape}')

	psf = psf / 128.0
	psf = np.roll(psf, U, axis=0)
	psf = np.roll(psf, U, axis=2)
	# print(psf)
	
	return psf

def resampling_operator(M):
	# Vars
	mtx = np.zeros((M**2, M))
	x = np.asarray([i for i in range(M**2)])
	xs = np.floor(x**0.5).astype(int)
	print(f'xs zeros: {len([i for i in xs if i == 0])}')
	mtx[x, xs] = 1

	xs[xs == 0] = 1
	mtx = scipy.sparse.csr_matrix(mtx)
	mtx = scipy.sparse.spdiags(1 / xs, 0, M**2, M**2) @ mtx 
	mtxi = mtx.T
	print(f'mtx shape: {mtx.shape}')
	print(f'mtxi shape: {mtxi.shape}')

	K = np.round(np.log(M) / np.log(2)).astype(int)
	print(f'K: {K}')
	for i in range(K):
		mtx = 0.5 * (mtx[::2, :] + mtx[1::2, :])
		mtxi = 0.5 * (mtxi[:, ::2] + mtxi[:, 1::2])

	return mtx.toarray(), mtxi.toarray()

def cnlos_reconstruction(mat_in_path):
	##### SETUP DATA #####
	# Constants
	isbackprop = False
	isdiffuse = False
	bin_resolution = 4e-12
	c = 3e8
	K = 2
	snr = 8e-1
	z_trim = 600
	z_offset_dict = {
		'data_resolution_chart_40cm.mat' : [350, 'psf.mat'],
		'data_resolution_chart_65cm.mat' : [700, 'psf.mat'],
		'data_dot_chart_40cm.mat' : [350, 'psf.mat'],
		'data_dot_chart_65cm.mat' : [700, 'psf.mat'],
		'data_mannequin.mat' : [300, 'psf.mat'],
		'data_exit_sign.mat' : [600, 'psf.mat'],
		'data_s_u.mat' : [800, 'psf.mat'],
		'data_outdoor_s.mat' : [700, 'psf_large.mat'],
		'data_diffuse_s.mat' : [100, 1, snr * 1e-1, 'psf.mat']
	}

	# Open matlab data file
	mat = scipy.io.loadmat(mat_in_path)
	mat_width, mat_rect_data = mat['width'][0, 0], mat['rect_data']
	print(f'mat_width: {mat_width}')
	print(f'mat_rect_data shape: {mat_rect_data.shape}')

	# Get fixed z offset and other constant params
	mat_fname = os.path.split(mat_in_path)[-1]
	print(f'mat_fname: {mat_fname}')
	mat_params = z_offset_dict[mat_fname]
	if len(mat_params) > 2: z_offset, isdiffuse, snr = mat_params[:3]
	else: z_offset = mat_params[0]
	psf_fname = mat_params[-1]
	print(f'z_offset, isdiffuse, snr: {z_offset}, {isdiffuse}, {snr}')

	# Get spatial, temporal, and range dims
	N, M = mat_rect_data.shape[0], mat_rect_data.shape[2]
	mat_range = M * c * bin_resolution
	print(f'M, N: {M}, {N}')

	# Downsample data to 16 picoseconds
	for i in range(K):
		M //= 2
		bin_resolution *= 2
		mat_rect_data = mat_rect_data[:, :, 0::2] + mat_rect_data[:, :, 1::2]
		z_trim = round(z_trim / 2)
		z_offset = round(z_offset / 2)
	print(f'M, N: {M}, {N}')

	# Set first group of histogram bins to zero (to remove direct component)
	mat_rect_data[:, : , 0:z_trim] = 0

	# Define NLOS blur kernel 
	psf = scipy.io.loadmat(psf_fname)['psf']
	print(f'psf shape: {psf.shape}')

	# Compute inverse filter of NLOS blur kernel
	fpsf = np.fft.fftn(psf)
	if isbackprop: invpsf = np.conj(fpsf)
	else: invpsf = np.conj(fpsf) / (np.abs(fpsf)**2 + 1.0 / snr)

	# Define volume representing voxel distance from wall
	Ml = np.linspace(0, 1, M).T
	grid_z = np.tile(Ml[:, None, None], [1, N, N])
	print(f'grid_z shape: {grid_z.shape}')

	# Get transform operators and permute data dims
	mtx, mtxi = resampling_operator(M)
	data = np.asarray(mat_rect_data.transpose(2, 1, 0))


	##### RUN ALGO #####
	# Step 1: Scale radiometric component
	if isdiffuse: data = data * grid_z**4
	else: data = data * grid_z**2
	print(f'data shape: {data.shape}')

	# Step 2: Resample time axis and pad result
	data_rs = data.reshape((data.shape[0], -1))
	print(f'data_rs shape: {data_rs.shape}')
	tdata = np.zeros((2 * M, 2 * N, 2 * N))
	tdata[:M, :N, :N] = (mtx @ data_rs).reshape((M, N, N))
	print(f'tdata shape: {tdata.shape}')


	# Step 3: Convolve with inverse filter and unpad result
	tvol = np.fft.ifftn(np.fft.fftn(tdata) * invpsf)
	tvol = tvol[:M, :N, :N]
	print(f'tvol shape: {tvol.shape}')

	# Step 4: Resample depth axis and clamp results
	tvol_rs = tvol.reshape((tvol.shape[0], -1))
	vol = (mtxi @ tvol_rs).reshape((M, N, N))
	vol = np.real(vol)
	vol[vol < 0] = 0
	print(f'vol shape: {vol.shape}')
	
	# Crop and flip for visualization
	ind = round(M * 2 * mat_width / (mat_range / 2))
	print(f'ind: {ind}, {M}, {mat_width}, {mat_range}, {z_offset}')
	vol = vol[:, :, -1::-1]
	vol = vol[z_offset:z_offset+ind, :, :]
	print(f'vol shape: {vol.shape}')
	vol1 = vol.max(0)
	print(f'vol1 shape: {vol1.shape}')

	# Display
	plt.figure('bob')
	plt.subplot(1, 3, 1)
	plt.imshow(vol1, cmap='gray')
	plt.show()

def main():
	# Args
	try:
		mat_in_path = sys.argv[1]
	except:
		print(f'\n***ERROR*** Must have positional argument for path to matlab file:\n\n$ python3 project.py path/to/data.m\n')
		sys.exit()

	# CNLOS Reconstruction
	cnlos_reconstruction(mat_in_path)

if __name__ == '__main__':
	main()