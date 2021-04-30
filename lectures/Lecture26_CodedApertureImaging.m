%% Basic coded-aperture imaging
clear;
clc;
test_im = double(rgb2gray(imread('./test_image_sq.png')))./255;
 
 
 
mask_code = zeros(size(test_im));
numPinholes = 6;
 
 
pinhole_locations = randi(prod(size(test_im)), numPinholes);
mask_code(pinhole_locations) = 1./numPinholes;
 
figure(1); clf;
subplot(221)
imagesc(mask_code); 
 
 
subplot(222)
imshow(test_im)
H_sys = fft2(mask_code);
 
measured_image = ifft2(H_sys.*fft2(test_im));
 
subplot(223)
imshow(measured_image./max(measured_image(:)))
 
 
subplot(224)
reconstruction = ifft2(fft2(measured_image)./H_sys);
imshow(reconstruction)
 
 
 
%% With noise
 
measured_image_noisy = awgn(measured_image,30, 'measured');
figure(2); clf;
subplot(221)
imshow(measured_image_noisy./max(measured_image_noisy(:)))
 
subplot(222)
reconstruction_noisy = ifft2(fft2(measured_image_noisy)./H_sys);
imshow(reconstruction_noisy)
 
lambda = .00005;
subplot(223)
tik_reconstruction_noisy = ifft2(((H_sys')./(H_sys' .* H_sys + lambda)).*fft2(measured_image_noisy) );
imshow(tik_reconstruction_noisy./max(tik_reconstruction_noisy));
 

