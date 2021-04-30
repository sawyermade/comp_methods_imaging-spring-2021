%% Lecture 9 - Matlab practice 6: FFT
%
%
% Example 1: consider signal:
%        f(t) = 0.7 sin (100*pi*t) + 2*sin(240*pi*t)
% (a) Sample the signal: obtain 1500 samples at a 
%     sampling frequency of 1000 Hz
% (b) Plot the magnitude spectrum
%

Fs = 1000;
T = 1/Fs;
L = 1500;

t = (0:L-1)*T;

signal_ft = 0.7*sin(100*pi*t) + 2*sin(240*pi*t);

figure(1);
plot(t, signal_ft)

figure(2);
signal_FW = fft(signal_ft)
signal_FW(100)
freq_axis = Fs*(-(L/2):(L/2 - 1))/L
plot(freq_axis, abs(fftshift(signal_FW))/L)




%%  Example 2: Displaying image and computing FT of image
%
% Load 'rocky.png' and 'rocky_brick.png' image;
% Display them;
% Compute FFT and display, comment on differences

rocky = imread('rocky.png');

rocky_brick = imread('rocky_brick.png');

figure(1)
subplot(121)
imshow(rocky)

subplot(122)
imshow(rocky_brick)

figure(2)
R1 = fft2(rocky(:,:,1));
RB = fft2(rocky_brick(:,:,1));

abs_R1 = fftshift(abs(R1));
abs_RB = fftshift(abs(RB));
%%
subplot(121)
imshow(abs_R1, [0 20000]);

subplot(122)
imshow(abs_RB, [0 20000])



















