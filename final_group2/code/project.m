% Default scene: SET TO 1 through 9
scene = 8;

% Constants
bin_resolution = 4e-12; % Native bin resolution for SPAD is 4 ps
c              = 3e8;   % Speed of light (meters per second)

% Adjustable parameters
isbackprop = 0;         % Toggle backprojection
isdiffuse  = 0;         % Toggle diffuse reflection
K          = 2;         % Downsample data to (4 ps) * 2^K = 16 ps for K = 2
snr        = 8e-1;      % SNR value
z_trim     = 600;       % Set first 600 bins to zero


% Load scene & set visualization parameter
switch scene
case {1}
load data/data_resolution_chart_40cm.mat
z_offset = 350;
case {2}
load data/data_resolution_chart_65cm.mat
z_offset = 700;
case {3}
load data/data_dot_chart_40cm.mat
z_offset = 350;
case {4}
load data/data_dot_chart_65cm.mat
z_offset = 700;
case {5}
load data/data_mannequin.mat
z_offset = 300;
case {6}
load data/data_exit_sign.mat
z_offset = 600;
case {7}
load data/data_s_u.mat
z_offset = 800;
case {8}
load data/data_outdoor_s.mat
z_offset = 700;
case {9}
load data/data_diffuse_s.mat
z_offset = 100;

% Because the scene is diffuse, toggle the diffuse flag and 
% adjust SNR value correspondingly.
isdiffuse = 1;
snr = snr.*1e-1;
end

N = size(rect_data,1);        % Spatial resolution of data
M = size(rect_data,3);        % Temporal resolution of data
range = M.*c.*bin_resolution; % Maximum range for histogram

% Downsample data to 16 picoseconds
for k = 1:K
    M = M./2;
    bin_resolution = 2*bin_resolution;
    rect_data = rect_data(:,:,1:2:end) + rect_data(:,:,2:2:end);
    z_trim = round(z_trim./2);
    z_offset = round(z_offset./2);
end

% Set first group of histogram bins to zero (to remove direct component)
rect_data(:,:,1:z_trim) = 0;

% Define NLOS blur kernel 
psf = definePsf(N,M,width./range);
% display('psf:')

% Compute inverse filter of NLOS blur kernel
fpsf = fftn(psf);
if (~isbackprop)
invpsf = conj(fpsf) ./ (abs(fpsf).^2 + 1./snr);
else
invpsf = conj(fpsf);
end

% Define transform operators
[mtx,mtxi] = resamplingOperator(M);

% Permute data dimensions
data = permute(rect_data,[3 2 1]);
data_og = data;

% Define volume representing voxel distance from wall
Ml = linspace(0,1,M)';
grid_z = repmat(Ml,[1 N N]);

display('Inverting...');
tic;

% Step 1: Scale radiometric component
if (isdiffuse)
data = data.*(grid_z.^4);
else
data = data.*(grid_z.^2);
end

% Step 2: Resample time axis and pad result
data_reshape = data(:,:);
tdata = zeros(2.*M,2.*N,2.*N);
tdata(1:end./2,1:end./2,1:end./2)  = reshape(mtx*data(:,:),[M N N]);

% Step 3: Convolve with inverse filter and unpad result
tvol = ifftn(fftn(tdata).*invpsf);
tvol = tvol(1:end./2,1:end./2,1:end./2);

% Step 4: Resample depth axis and clamp results
vol  = reshape(mtxi*tvol(:,:),[M N N]);
vol_og = vol;
vol  = max(real(vol),0);
vol_0 = vol;


display('... done.');
time_elapsed = toc;

display(sprintf(['Reconstructed volume of size %d x %d x %d '...
'in %f seconds'], size(vol,3),size(vol,2),size(vol,1),time_elapsed));

tic_z = linspace(0,range./2,size(vol,1));
tic_y = linspace(-width,width,size(vol,2));
tic_x = linspace(-width,width,size(vol,3));

% Crop and flip reconstructed volume for visualization
ind = round(M.*2.*width./(range./2));
vol = vol(:,:,end:-1:1);
vol = vol((1:ind)+z_offset,:,:);
vol1 = max(vol, [], 1);

tic_z = tic_z((1:ind)+z_offset);

% View result
figure('pos',[10 10 900 300]);

subplot(1,3,1);
imagesc(tic_x,tic_y,squeeze(max(vol,[],1)));
title('Front view');
set(gca,'XTick',linspace(min(tic_x),max(tic_x),3));
set(gca,'YTick',linspace(min(tic_y),max(tic_y),3));
xlabel('x (m)');
ylabel('y (m)');
colormap('gray');
axis square;

subplot(1,3,2);
imagesc(tic_x,tic_z,squeeze(max(vol,[],2)));
% imagesc(tic_x,tic_z,vol1_mine);
title('Top view');
set(gca,'XTick',linspace(min(tic_x),max(tic_x),3));
set(gca,'YTick',linspace(min(tic_z),max(tic_z),3));
xlabel('x (m)');
ylabel('z (m)');
colormap('gray');
axis square;

subplot(1,3,3);
imagesc(tic_z,tic_y,squeeze(max(vol,[],3))')
title('Side view');
set(gca,'XTick',linspace(min(tic_z),max(tic_z),3));
set(gca,'YTick',linspace(min(tic_y),max(tic_y),3));
xlabel('z (m)');
ylabel('y (m)');
colormap('gray');
axis square;

% save('vars_mannequin.mat')

function psf = definePsf(U,V,slope)
    % Local function to computeD NLOS blur kernel

    x = linspace(-1,1,2.*U);
    y = linspace(-1,1,2.*U);
    z = linspace(0,2,2.*V);
    [grid_z,grid_y,grid_x] = ndgrid(z,y,x);

    % Define PSF
    psf = abs(((4.*slope).^2).*(grid_x.^2 + grid_y.^2) - grid_z);
%     psf = min(psf, [], 1)
    psf = double(psf == repmat(min(psf,[],1),[2.*V 1 1]));
    psf_sum = sum(psf(:, U, U));
    display('psf_sum: ')
    psf_sum
    psf = psf./sum(psf(:,U,U));
    psf_norm = norm(psf(:));
    psf_norm
    psf = psf./norm(psf(:));
    psfog = psf;
    psf = circshift(psf,[0 U U]);
end

function [mtx,mtxi] = resamplingOperator(M)
    % Local function that defines resampling operators

    mtx = sparse([],[],[],M.^2,M,M.^2);

    x = 1:M.^2;
    sx = sqrt(x);
    mtx_size = size(mtx);
    mtx_sub = sub2ind(size(mtx), x, ceil(sqrt(x)));
    mtx(sub2ind(size(mtx),x,ceil(sqrt(x)))) = 1;
    mtx_sum = sum(mtx, 'all'); 
    mtx_diag = spdiags(1./sqrt(x)',0,M.^2,M.^2);
    mtx  = spdiags(1./sqrt(x)',0,M.^2,M.^2)*mtx;
    mtxi = mtx';

    K = log(M)./log(2);
    for k = 1:round(K)
        mtx  = 0.5.*(mtx(1:2:end,:)  + mtx(2:2:end,:));
        mtxi = 0.5.*(mtxi(:,1:2:end) + mtxi(:,2:2:end));
    end
end
