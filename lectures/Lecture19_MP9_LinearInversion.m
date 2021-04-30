%% Lecture 19: Matlab practice 9 - Linear Inversion

%% Part I: Deconvolution/Deblurring Problem - Forward model

clear all;
clc;
 
% Important functions
FT_2d = @(x) fftshift(fft2(ifftshift(x)));
IFT_2d = @(x) fftshift(ifft2(ifftshift(x)));
 
 
N = 1000;
M = 1000;
 
Delta = 2;
delta = 2;
 
sigma2 = 5000;
 
dw = 1/N/Delta;
omega_x = (-(N/2):(N/2 -1))*dw;
omega_y = (-(N/2):(N/2 -1))*dw;
 
[omega_xx, omega_yy ] = meshgrid(omega_x, omega_y);
 
% Fourier transform of the object basis
P_00 = sinc(Delta*omega_xx).*sinc(Delta*omega_yy);
 
% Pixel sampling transfer function
PSI_00 = P_00;
 
% Compute TF of the optical part (OTF)
H = 2*pi*sigma2*exp( -0.5*sigma2*(omega_xx.^2 + omega_yy.^2) );
 
% Compute overall system transfer function
H_sys = P_00 .* PSI_00 .* H;

figure(1)
imagesc(omega_x, omega_y, real(H_sys))

%%
figure(2);
subplot(121)
ff_in = imread('test_image_sq.png');
ff_in = rgb2gray(ff_in);
imshow(ff_in)
title('Original')

subplot(122)
FF_out_blurred = FT_2d(ff_in).* H_sys;
ff_out_blurred = IFT_2d(FF_out_blurred);
ff_out_blurred = ff_out_blurred./max(ff_out_blurred(:));

imshow(ff_out_blurred)
title('System output (blurred)')

%% Additive white gaussian noise (AWGN)
SNR = 0;   %Signal to noise ratio
ff_out_blurred_noisy = awgn(ff_out_blurred, SNR, 'measured');

figure(3) 
imshow(ff_out_blurred_noisy)
title('System output (blurred + noise)')


%% Part II: Linear Inversion
%% (a): Unregularized - Standard Inverse and Pseudo-Inverse
%
% gout = h convolved with gin <---> Gout = H.*Gin
FF_out_bn = FT_2d(ff_out_blurred_noisy);
FF_in_est  = FF_out_bn./H_sys;
ff_in_est = IFT_2d(FF_in_est);
figure(4)
imshow(ff_in_est./max(ff_in_est(:)))

%% (b): Regularized inversion - Tikhonov
%

lambda = 0.00001;

FF_est_tik = ((H_sys')./(H_sys' .* H_sys + lambda)) .* FF_out_bn;
ff_est_tik = IFT_2d(FF_est_tik);

figure(5)
imshow(ff_est_tik./max(ff_est_tik(:)))
