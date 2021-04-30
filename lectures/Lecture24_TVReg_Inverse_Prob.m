%% TV regularized inversion of (discrete) linear forward models (g = Af + n)

% [reconstruction] = fistatvrecon(A, measurement, lambda_val, numpixels_xy);

% min_f ||Af - g||_2^2 + lambda_val * TV(f)

clear;
close all;
clc;

% Important functions
FFT_2d = @(x) fftshift(fft2(ifftshift(x)));
IFFT_2d = @(x) fftshift(ifft2(ifftshift(x)));

ff = double(imread('test_image_300x300.png'))/255;
ff = rgb2gray(ff);

ai = [1 4 6 9 14 20 32 44 60 44 32 20 14 9 6 4 1]';
ai = ai*ai';
ai = ai./sum(ai(:));

%
AA = convmtx2(ai, size(ff));
gg = AA*ff(:);

gg = reshape(gg,size(ff)+size(ai)-1);
%% TV regularized inversions


SNR = 20;
meas = awgn(gg, SNR, 'measured');
lambda_val = .003;

numpixels_xy = size(ff);

tic
ff_est_tv = fistatvrecon(AA, meas, lambda_val, numpixels_xy);
toc

%% Tikhonov
% Using fft techniques to compute the Tikhonov-reg inverse
tik_reg = 1e-1;
tic
FTA = FFT_2d(padarray(ai,[150 150]));
FTA = FTA(1:end-1,1:end-1);
 
GG = FFT_2d(gg);
FT_meas = FFT_2d(meas);
f_est_tik = IFFT_2d( ((FTA').* FT_meas)./(abs(FTA).^2 + tik_reg) );
toc


1%%
figure(1)
imshow(ff)
title('Ground truth')
figure(2)
imshow(gg)
title('Measurement')
figure(3)
imshow(ff_est_tv)
title('TV Reconstruction')
figure(4)
imshow(f_est_tik)
title('Tikhonov reconstruction')