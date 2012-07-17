% function []=test_grm_creation(nx,dx,nz,dz, a)
%% This m-file is a test script to explore the difference between different ways of generating a gaussain random media. 
%
%
%% setup default values for testing
% if nargin ==0
    nx = 1000;          % number of nodes in x direction
    dx = 20;            % x-spacing
    nz = 1000;          % number of nodes in z direction
    dz = 20;            % z-spacing
    a = 500;            % length of coprrelation filter
    std = 1500          % standard deviation of medium
    mean_medium = 5000; % mean of final medium
% end



%% Test 1
% 1) generate a grid x and assign a uniform random number between 0 and 1 to each grid point
% 2) take the fft(x)
% 3) convolve with the power spectra of the auto-cortrelation function I*x
% 4) take the ifft
% 5) scale with the desired std dev. (perturbation)
% 6) add the mean

randn('state',0)
rand('state',0)
% compute the positive wave numbers. 
kx = repmat([0:nx/2]./(nx*dx),[nz/2+1,1]);
kz = repmat([0:nz/2]'./(nz*dz),[1,nx/2+1]);
k = sqrt(kx.^2+kz.^2);

figure
med1 = rand(nz,nx);
subplot(4,2,1)
imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz] ,med1)
title('med1: uniform random numbers on grid')
colormap bone
subplot(4,2,2)
hist(med1(:))
title(['hist med1: \mu = ', num2str(mean(med1(:))), ' \sigma = ', num2str(sqrt(var(med1(:))))])

fft_med1 = fft2(med1);
fG = exp(-a^2.*k.^2/8);   % Power spectra of filter (+ve wave numbers)
fG = extend2negwave(fG,[nx,1,nz]); % Extend filter into (-ve wavenumbers)
subplot(4,2,3)
h= imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz], fG)
title('filter: ')
subplot(4,2,4)
colorbar(h, 'WestOutside')
axis off

fft_med2 = fft_med1.*fG; %fft of the desired medium
med2 = ifft2(fft_med2); %Desired medium
subplot(4,2,5)
imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz],real(med2));
colormap bone
title('med2')
subplot(4,2,6)
hist(real(med2(:)))
title(['med2 : \mu = ', num2str(mean(med2(:))), ' \sigma = ', num2str(sqrt(var(med2(:))))])

% let's look for imaginary values in med2
max_imag = max(abs(imag(med2(:))))
min_imag = min(abs(imag(med2(:))))

%Finally we scale for std and add the mean
med3 = med2 + mean_medium;
subplot(4,2,7)
imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz],real(med3))
title(['med3 = med2+mean '])
colormap(bone)
subplot(4,2,8)
hist(med3(:))
title(['med 3: \mu = ', num2str(mean(med3(:))), ' \sigma = ', num2str(sqrt(var(med3(:))))])  

test_med1 = med1;
test_fG = fG;
test_med2 = med2;
test_med3 = med3;

save test.mat test_med1 test_fG test_med2 test_med3


%% Test 2
% 1) generate a grid x and assign a normall distributed random number with mean=0 and std=1 to each grid point
% 2) take the fft(x)
% 3) convolve with the power spectra of the auto-cortrelation function I*x
% 4) take the ifft
% 5) scale with the desired std dev. (perturbation)
% 6) add the mean

randn('state',0)
rand('state',0)
% compute the positive wave numbers. 
kx = repmat([0:nx/2]./(nx*dx),nz/2+1,1);
kz = repmat([0:nz/2]'./(nz*dz),1,nx/2+1);
k = sqrt(kx.^2+kz.^2);

figure
med1 = randn(nz,nx);
subplot(4,2,1)
imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz] ,med1)
title('med1: gaussian random numbers on grid')
colormap bone
subplot(4,2,2)
hist(med1(:))
title('hist med1')

fft_med1 = fft2(med1);
fG = exp(-a^2.*k.^2/8);   % Power spectra of filter (+ve wave numbers)
fG = extend2negwave(fG,[nx,1,nz]); % Extend filter into (-ve wavenumbers)
subplot(4,2,3)
imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz], fG)
title('filter: ')
subplot(4,2,4)
colorbar(h, 'WestOutside')
axis off

fft_med2 = fft_med1.*fG; %fft of the desired medium
med2 = ifft2(fft_med2); %Desired medium
subplot(4,2,5)
imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz],real(med2));
colormap bone
title('med2')
subplot(4,2,6)
hist(real(med2(:)))
title('med2')

% let's look for imaginary values in med2
max_imag = max(abs(imag(med2(:))))
min_imag = min(abs(imag(med2(:))))

%Finally we scale for std and add the mean
med3 = med2 + mean_medium;
subplot(4,2,7)
imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz],real(med3))
title('med3 = med2+mean')
colormap(bone)
subplot(4,2,8)
hist(med3(:))
title('med 3')  



%% Test 3
% 1) generate a grid x and assign a uniform random number r to each grid value: r in (0,1)
% 2) compute a random phase at each grip point using exp(i2*pi*r)  
% 3) convolve with the power spectra of the auto-cortrelation function I*x (note that I includes the fractional varience epsilon)
% 4) take the ifft
% 5) add the mean
clear kx kz k fft_med1 fft_med2 med1 med2 med3
randn('state',0)
rand('state',0)
% compute the positive wave numbers. 
kx = repmat([0:nx/2]./(nx*dx),nz/2+1,1);
kz = repmat([0:nz/2]'./(nz*dz),1,nx/2+1);
k = sqrt(kx.^2+kz.^2);

figure
clear i
med1 = exp(i*pi*2*rand(nz,nx));
%subplot(4,2,1)
% imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz] ,med1)
% title('med1: uniform random numbers on grid')
% colormap bone
% subplot(4,2,2)
% hist(med1(:))
% title('hist med1 -exp(2*pi*rand(nx,nz))')


fG = exp(-a^2.*k.^2/8);   % Power spectra of filter (+ve wave numbers)
fG = extend2negwave(fG,[nx,1,nz]); % Extend filter into (-ve wavenumbers)
subplot(4,2,3)
imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz], fG)
title('filter: ')
subplot(4,2,4)
colorbar(h, 'WestOutside')
axis off

fft_med2 = med1.*fG; %fft of the desired medium (note no need to fft(med1))
med2 = ifft2(fft_med2); %Desired medium
subplot(4,2,5)
imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz],real(med2));
colormap bone
title('med2')
subplot(4,2,6)
hist(real(med2(:)))
title('med2')

% let's look for imaginary values in med2
max_imag = max(abs(imag(med2(:))))
min_imag = min(abs(imag(med2(:))))

%Finally we scale for std and add the mean
med3 = med2 + mean_medium;
subplot(4,2,7)
imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz],real(med3))
title('med3 = med2+mean')
colormap(bone)
subplot(4,2,8)
hist(med3(:))
title('med 3')  

figure
subplot(2,2,1)
hist(real(med2(:)))
title('hist(real(med2(:)))')
subplot(2,2,2)
hist(imag(med2(:)))
title('hist(imag(med2(:)))')
subplot(2,2,3)
hist(real(med3(:)))
title('hist(real(med3(:)))')
subplot(2,2,4)
hist(imag(med3(:)))
title('hist(imag(med3(:)))')


%%======================




