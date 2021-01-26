% Beam Simulation Demo (Version 14.0)
%
% %%%%%%%%%%%%%%%%%%%%%         Input        %%%%%%%%%%%%%%%%%%
% 1. Volume Size (comprise several planes)
% 2. Beam Type; Beam Parameter; 
% 3. Function Options
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%% Packages & Functions %%%%%%%%%%%%%%%%%%
% 0. Complex amplitude of boundaries
% 1. Three dimensional complex amplitude distributions (Two light beams)
% 2. Complex amplitude distributions of a certain plane
% 3. Average Poynting vectors of a certain plane
% 4. Average angular momentum distributions of a certain plane
% 5. Interference pattern of a certain plane
% 6. DMD holograms (Two light beams)
% 7. Experimental simulations of diffractions
% 8. Isosurface of intensity
% 9. 3D GS Algorithm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by: Yue chen
% Department of Modern Physics, University of Science and Technology of
% China (USTC)
% December. 4, 2018
% More imformation£ºchenyue8@mail.ustc.edu.cn
% Before simulation, one need to check readme.txt.

% parameters: Sub-nyquist sampling principle; zR; Paths;

clc ; clear all ; close all ; dbstop if error ;

%%  0. Input Setup
% 0. Add path
path(path,genpath(pwd));
% 1.Global invariants
% a,b,c for size of volume; M,N,MN for samples of volume;(SI UNINS)
global a ; global b ; global c ; global M ; global N ; global MN ;
global xa ; global ya ; global za ; global x ; global y ; global z ; global r ;
a = 1e-2 ; M = 1080 ; b= 1e-2 ; N = 1080 ; c = 1.25e-3 ; MN = 3 ;
xa = linspace ( -a , a , M ) ; ya = linspace ( b , -b , N ) ; za = linspace ( -c , c , MN ) ;
[ x , y , z ] = meshgrid ( xa , ya , za ) ; r = sqrt ( x.^2 + y.^2 ) ;
% 2.Input beams
[ C0 , Y0 ] = GPOV () ;

%%  1. Three dimensional complex amplitude distributions (Two light beams)
% 1. Beam Propagation Method ( Options: 0 for non-focusing diffraction & 1
% for focusing diffraction; 0 for D-FFT & 1 for S-FFT) 
[ C , Czi , Cxoz ] = Diffraction ( C0 , 0 ) ;
[ Y , Yzi , Yxoz ] = Diffraction ( Y0 , 0 ) ;
% 2. Figure 1
Fig1 ( C0 , Czi , Cxoz , Y0 , Yzi , Yxoz ) ;

%  2. DMD Logistic Holograms (Two light beams)
% 1.SuperPixel Method (Reference: ) 
[ DMDpixels_C , DMDpixels_Y ] = SuperPixelMethod ( C0 , Y0 , [] ) ;
% % 3. Binary Mask
% DMDpixels1_logical = logical ( abs ( DMDpixels_C ) ) ;
% DMDpixels2_logical = logical ( abs ( DMDpixels_Y ) ) ;
% 4. Figure
Fig2 ( DMDpixels_C , DMDpixels_Y , [] ) ; % gray
% Fig2 ( DMDpixels1_logical , DMDpixels2_logical , [] ) ; % logical






