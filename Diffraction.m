%  Beam Propagation in free space
% 
% %%%%%%%%%%%%%%%%%%%%%         Description        %%%%%%%%%%%%%%%%%%
% 1. Calculate the three dimensional optical field via D-FFT method.
% 2. The wavelength and focoal length should also be changed while applying
%    in a practical optical configuration. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by: Yue chen
% Department of Modern Physics, University of Science and Technology of
% China (USTC)
% December. 4, 2018
% More imformation£ºchenyue8@mail.ustc.edu.cn

function [ E , Ezi , Exoz ] = Diffraction ( E0 )
% 1. Environrment setup
% Global invariants
global a;global b;global c;global M;global N;global MN;
global xa;global ya;global za;global x;global y;global z;global r;
lamda = 632.8e-9 ;
% 2. Diffraction 
E = zeros ( N , M , MN ) ;
% G0=fftshift(fft2(ifftshift(E0)));
G0 = fftshift ( fft2 ( E0 ) ) ;
% Initial plane excepted
k0 = 1 + floor ( MN / 2 ) ;
E ( : , : , k0 ) = E0 ;
% Beam propagation method
% Negetive Position calculation
for k = 1 : floor ( MN / 2 )
    % D-FFT method (caution for opposite calculation due to the huge value)
    Gz_conj = conj(G0) .* exp ( 1i * 2 * pi * abs ( z ( 1 , 1 , k ) ) .*...
        sqrt ( 1 / ( lamda .^ 2 ) -...
        ( x ( : , : , 1 ) / a / a / 4 * ( M - 1 ) ) .^ 2 - ( y ( : , : , 1 )...
        / b / b / 4 * ( N - 1 ) ) .^ 2 ) ) ;
%     Gz=G0.*exp(1i*2*pi*z(1,1,k)/lamda.*(1-lamda^2/2*((x/4/a/a*M).^2+(y/4/b/b*N).^2)));
%     Ez=fftshift(ifft2(ifftshift(Gz)));
    % Circ Filter
    Gz = conj ( Gz_conj ) ;
    % When using indicators, circle is fine; Using filter, ellipse.
    % Not necessary due to the evanescent wave 
    T = CircFilter ( 1 / lamda * 4 * a * a / ( M - 1 ) , 1 / lamda * 4 * b...
        * b / ( N - 1 ) ) ; 
    Gz = Gz .* T ;
    Ez = ifft2 ( ifftshift ( Gz ) ) ;
    % Normalized ( Generally unnecessary)
    % Ez=Ez./max(abs(Ez(:)));
    for i = 1 : N
        for j = 1 : M
% Fill Ez into Ezz
          E ( i , j , k ) = Ez ( i , j ) ;     
        end
    end
end
% Positive calculation
for k = ( floor ( MN / 2 ) + 2 ) : MN
    % D-FFT method (caution for opposite calculation due to the huge value)
    Gz = G0 .* exp ( 1i * 2 * pi * abs ( z ( 1 , 1 , k ) ) .* sqrt ( 1 /...
        ( lamda .^ 2 ) -...
        ( x ( : , : , 1 ) / a / a / 4 * ( M - 1 ) ) .^ 2 - ( y ( : , : , 1 )...
        / b / b / 4 * ( N - 1 ) ) .^ 2 ) ) ;
%     Gz=G0.*exp(1i*2*pi*z(1,1,k)/lamda.*(1-lamda^2/2*((x/4/a/a*M).^2+(y/4/b/b*N).^2)));
%     Ez=fftshift(ifft2(ifftshift(Gz)));
    % Circ Filter
    T = CircFilter ( 1 / lamda * 4 * a * a / ( M - 1 ) , 1 / lamda * 4 * b...
        * b / ( N - 1 ) ) ;
    Gz = Gz .* T ;
    Ez = ifft2 ( ifftshift ( Gz ) ) ; 
    % Normalized
    % Ez=Ez./max(abs(Ez(:)));
    for i = 1 : N
        for j = 1 : M
% Fill Ez into Ezz
          E ( i , j , k ) = Ez ( i , j ) ;     
        end
    end
end
% 3.Ezi
Ezi = E ( : , : , MN ) ;
% 4.xoz 
Dxoz = floor ( N / 2 ) ;
Exoz = zeros ( M , MN ) ;
for j = 1 : M
    for k = 1 : MN
        Exoz ( j , k ) = E (Dxoz , j , k ) ;
    end
end
end