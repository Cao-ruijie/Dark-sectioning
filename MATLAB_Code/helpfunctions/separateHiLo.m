function [Hi,Lo,lp,EL] = separateHiLo(image,params,deg,divide)

%% 基本参数
Nx = params.Nx;
Ny = params.Ny;
NA = params.NA;
emwavelength = params.emwavelength;
pixel_size = params.pixelsize;

%% 其他参数
res = 0.5 * emwavelength / NA/ params.factor;     % resolution
k_m = Ny / (res / pixel_size(1));    % objective cut-off frequency ???
kc = floor(k_m * 0.2);             % cut-off frequency between hp and lp filter
sigmaLP = kc*2/2.355;                % Finding sigma value for low pass

%% 滤波
lp = lpgauss(Nx,Ny,sigmaLP*2*divide);
hp = hpgauss(Nx,Ny,sigmaLP*2*divide);
elp = lpgauss(Nx,Ny,sigmaLP/deg);
ehp = hpgauss(Nx,Ny,sigmaLP/deg);

%% 得到高低频和极低频率
Hi = real(ifft2(fft2(image).*hp));
Lo = real(ifft2(fft2(image).*lp));

fft_image = fftshift(fft2(image));
Hi = real(ifft2(ifftshift(fft_image.*fftshift(hp))));
Lo = real(ifft2(ifftshift(fft_image.*fftshift(lp))));

% elp = zeros(Nx,Ny);
% elp(floor((Nx+1)/2)-1:floor((Nx+1)/2)+1,floor((Nx+1)/2)-1:floor((Nx+1)/2)+1) = 1;
EL = real(ifft2(fft2(image).*elp));
EH = real(ifft2(fft2(image).*ehp));
end

function [ out ] = hpgauss(H,W,SIGMA)
%   Creates a 2D Gaussian filter for a Fourier space image of height H and
%   width W. SIGMA is the standard deviation of the Gaussian.
out=1-lpgauss(H,W,SIGMA);
end

function [ out ] = lpgauss(H,W,SIGMA)
%   Creates a 2D Gaussian filter for a Fourier space image
%   W is the number of columns of the source image and H is the number of
%   rows. SIGMA is the standard deviation of the Gaussian.
H = double(H);
W = double(W);
kcx = (SIGMA);
kcy = ((H/W)*SIGMA);
temp0 = -floor(W/2);
[x,y] = meshgrid(-floor(W/2):floor((W-1)/2), -floor(H/2):floor((H-1)/2));
temp = -(x.^2/(kcx^2)+y.^2/(kcy^2));
out = ifftshift(exp(temp));
% out = ifftshift(exp(-(x.^2/(kcx^2)+y.^2/(kcy^2))));
end
