%  Superpixel Method for DMD (DMD: 1920*1080)
% (Reference: Sebastianus, A, Goorden, et al. 
%   Superpixel-based spatial amplitude and phase modulation using
%   a digital micromirror device[J]. Optics Express, 2014.) 

% %%%%%%%%%%%%%%%%%%%%%         Description        %%%%%%%%%%%%%%%%%%
% 1. Calculate the hologram for DMD via superpixel method
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by: Yue chen
% Department of Modern Physics, University of Science and Technology of
% China (USTC)
% December. 4, 2018
% More imformation£ºchenyue8@mail.ustc.edu.cn

function [DMDpixels1,DMDpixels2]=SuperPixelMethod(E1,E2,superpixelSize)
% 1.²ÎÊý
global a;global b;global c;global M;global N;global MN;
global xa;global ya;global za;global x;global y;global z;global r;
if isempty(superpixelSize)==1
    superpixelSize=4;
end
% 2.DMDMasks
[NN,MM]=size(E1);
maskCenterX = ceil((MM+1)/2);
maskCenterY = ceil((NN+1)/2);
FourierMaskSuperpixelResolution = zeros(NN,MM);
nSuperpixelX = MM / superpixelSize;
nSuperpixelY = NN / superpixelSize;
for i = 1:NN
    for j = 1:MM
        if (i-maskCenterY)^2+(NN/MM*(j-maskCenterX))^2 < (NN/superpixelSize/2)^2
            FourierMaskSuperpixelResolution(i,j) = 1;
        end
    end
end
for i=1:2
    if i==1
E_target_ft = fftshift(fft2(ifftshift(E1)));
% Apply mask
E_target_ft = FourierMaskSuperpixelResolution.*E_target_ft;
% Remove zeros outside of mask
E_superpixelResolution_ft = E_target_ft((maskCenterY - ceil((nSuperpixelY-1)/2)):(maskCenterY + floor((nSuperpixelY-1)/2)),(maskCenterX - ceil((nSuperpixelX-1)/2)):(maskCenterX + floor((nSuperpixelX-1)/2)));
% Add phase gradient to compensate for anomalous 1.5 pixel shift in real
% plane
phaseFactor = ones(nSuperpixelY,nSuperpixelX);
for k = 1:nSuperpixelY
    for j = 1:nSuperpixelX
        phaseFactor(k,j) = exp(2*1i*pi*(k/nSuperpixelY+j/nSuperpixelX)*3/8);
    end
end
E_superpixelResolution_ft = E_superpixelResolution_ft.*phaseFactor;
% Fourier transform back to DMD plane
E_superpixelResolution = fftshift(ifft2(ifftshift( E_superpixelResolution_ft )));
% Normalize such that the maximum amplitude is 1
E_superpixelResolution = E_superpixelResolution / max(max(abs(E_superpixelResolution)));
m = 4;
normalization='maxAmplitude';
% Load lookup table
targetFields = dlmread('targetFields');
lookupTable = dlmread('lookupTable');
gridParameters = dlmread('gridParameters');
maxAmplitude = gridParameters(1);
stepSize = gridParameters(2);
lookupTable_x0 = (length(lookupTable)+1)/2;
[ny,nx] = size(E_superpixelResolution);
DMDpixels = zeros(ny*m,nx*m);
% Decrease maximum amplitude to 1 if needed
E_superpixelResolution= exp(angle(E_superpixelResolution)*1i).*min(abs(E_superpixelResolution),1);
% Correct for overall phase offset
E_superpixelResolution = E_superpixelResolution.*exp(11/16*pi*1i);
% Choose normalization: maxAmplitude for highest efficiency, highRes to
% restrict the modulation to a smaller and denser disk in the complex plane
switch normalization
    case 'maxAmplitude'
        E_superpixelResolution =E_superpixelResolution * maxAmplitude;
    case 'highRes'
        E_superpixelResolution=E_superpixelResolution * 0.906131;       % 4 pixel, to be calibrated for different superpixel sizes
end
% Loop over superpixels. Find correct combination of pixels to turn on in the lookup table and put
% them into the 'DMDpixels' matrix that contains the DMD pattern
E_superpixelResolution= E_superpixelResolution / stepSize;
for j = 0:ny-1
    for i = 0:nx-1
        idx = lookupTable(round(imag(E_superpixelResolution(j+1,i+1))+lookupTable_x0),round(real(E_superpixelResolution(j+1,i+1))+lookupTable_x0));
        pixels = targetFields(idx,2:m^2+1);
        shift = mod(m*j,m^2);
        pixels = [pixels(shift+1:m^2) pixels(1:shift)];
        DMDpixels(m*j+1:m*j+m,m*i+1:m*i+m) = reshape(pixels,m,m);
    end
end
% Apply phase gradient to DMD pixels (equivalent to placing Fourier mask
% off-axis):
phaseFactor = ones(ny*m,nx*m);
for k = 1:ny*m
    for j = 1:nx*m
        phaseFactor(k,j) = exp(1i*pi*(k+4*j)/8);
    end
end
DMDpixels = DMDpixels.*phaseFactor;
DMDpixels1=DMDpixels;
    else
E_target_ft = fftshift(fft2(ifftshift(E2)));
% Apply mask
E_target_ft = FourierMaskSuperpixelResolution.*E_target_ft;
% Remove zeros outside of mask
E_superpixelResolution_ft = E_target_ft((maskCenterY - ceil((nSuperpixelY-1)/2)):(maskCenterY + floor((nSuperpixelY-1)/2)),(maskCenterX - ceil((nSuperpixelX-1)/2)):(maskCenterX + floor((nSuperpixelX-1)/2)));
% Add phase gradient to compensate for anomalous 1.5 pixel shift in real
% plane
phaseFactor = ones(nSuperpixelY,nSuperpixelX);
for k = 1:nSuperpixelY
    for j = 1:nSuperpixelX
        phaseFactor(k,j) = exp(2*1i*pi*(k/nSuperpixelY+j/nSuperpixelX)*3/8);
    end
end
E_superpixelResolution_ft = E_superpixelResolution_ft.*phaseFactor;
% Fourier transform back to DMD plane
E_superpixelResolution = fftshift(ifft2(ifftshift( E_superpixelResolution_ft )));
% Normalize such that the maximum amplitude is 1
E_superpixelResolution = E_superpixelResolution / max(max(abs(E_superpixelResolution)));
m = 4;
normalization='maxAmplitude';
% Load lookup table
targetFields = dlmread('targetFields');
lookupTable = dlmread('lookupTable');
gridParameters = dlmread('gridParameters');
maxAmplitude = gridParameters(1);
stepSize = gridParameters(2);
lookupTable_x0 = (length(lookupTable)+1)/2;
[ny,nx] = size(E_superpixelResolution);
DMDpixels = zeros(ny*m,nx*m);
% Decrease maximum amplitude to 1 if needed
E_superpixelResolution= exp(angle(E_superpixelResolution)*1i).*min(abs(E_superpixelResolution),1);
% Correct for overall phase offset
E_superpixelResolution = E_superpixelResolution.*exp(11/16*pi*1i);
% Choose normalization: maxAmplitude for highest efficiency, highRes to
% restrict the modulation to a smaller and denser disk in the complex plane
switch normalization
    case 'maxAmplitude'
        E_superpixelResolution =E_superpixelResolution * maxAmplitude;
    case 'highRes'
        E_superpixelResolution=E_superpixelResolution * 0.906131;       % 4 pixel, to be calibrated for different superpixel sizes
end
% Loop over superpixels. Find correct combination of pixels to turn on in the lookup table and put
% them into the 'DMDpixels' matrix that contains the DMD pattern
E_superpixelResolution= E_superpixelResolution / stepSize;
for j = 0:ny-1
    for i = 0:nx-1
        idx = lookupTable(round(imag(E_superpixelResolution(j+1,i+1))+lookupTable_x0),round(real(E_superpixelResolution(j+1,i+1))+lookupTable_x0));
        pixels = targetFields(idx,2:m^2+1);
        shift = mod(m*j,m^2);
        pixels = [pixels(shift+1:m^2) pixels(1:shift)];
        DMDpixels(m*j+1:m*j+m,m*i+1:m*i+m) = reshape(pixels,m,m);
    end
end
% Apply phase gradient to DMD pixels (equivalent to placing Fourier mask
% off-axis):
phaseFactor = ones(ny*m,nx*m);
for k = 1:ny*m
    for j = 1:nx*m
        phaseFactor(k,j) = exp(1i*pi*(k+4*j)/8);
    end
end
DMDpixels = DMDpixels.*phaseFactor;
DMDpixels2=DMDpixels;
    end
end
end