% Generalized perfect optical vortex beams
% 
% %%%%%%%%%%%%%%%%%%%%%         Description        %%%%%%%%%%%%%%%%%%
% 1. The expression of GPOVs
% 2. One may change the expression of curves and topololgical of beams if
%    needed.
% 3. The wavelength and focoal length should also be changed while applying
%    in a practical optical configuration. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by: Yue chen
% Department of Modern Physics, University of Science and Technology of
% China (USTC)
% December. 4, 2018
% More imformation£ºchenyue8@mail.ustc.edu.cn

function [ C , Y ] = GPOV()
% 1.
global a ; global b;global c;global M;global N;global MN;
global xa;global ya;global za;global x;global y;global z;global r;
% 2. 3D curves

% Ring
R0 = 1e-4 ;
T = 2 * pi ;
P = 1001 ;
t = linspace ( 0 , T , P ) ;
Curve_x = R0 .* cos ( t ) ; 
Curve_y = R0 .* sin ( t ) ;

% % Ellipse curve
% R0 = 0.5e-4 ;
% aa = R0 ;
% bb = 0.5 * aa ;
% T = 2 * pi ;
% P = 1001 ;
% t = linspace ( 0 , T , P ) ;
% Curve_x = aa .* cos ( t ) ;
% Curve_y = bb .* sin ( t ) ;

% % Astroid 1.5e-4
% R0 = 2e-4 ;
% T = 2*pi ;
% P = 1001 ;
% t = linspace ( 0 , T , P ) ;
% Curve_x = 3 * R0 .* cos ( t ) + R0 .* cos ( 3 * t ) ;
% Curve_y = 3 * R0 .* sin ( t ) - R0 .* sin ( 3 * t ) ;

% % Archimde spiral
% R0 = 5.6e-4 ;
% T = 1 ;
% P = 1001 ;
% t = linspace ( 0 , T , P ) ;
% Curve_x = R0 .* t .* cos ( 10 * t ) ;
% Curve_y = R0 .* t .* sin ( 10 * t ) ;

% % creatative curve from mathematica
% R0=5.6e-4;
% T=2*pi;
% P=1001;
% t=linspace(0,T,P);
% Curve_x=0.27132E3+(-0.427983E2).*cos(t)+0.188097E2.*cos(2.*t)+( ...
%   -0.184248E2).*cos(3.*t)+0.107918E2.*cos(4.*t)+0.346443E1.*cos(5.* ...
%   t)+0.153664E1.*cos(6.*t)+0.746E1.*cos(7.*t)+(-0.287869E1).*cos(8.* ...
%   t)+(-0.864026E1).*cos(9.*t)+(-0.262604E1).*cos(10.*t)+0.361057E1.* ...
%   cos(11.*t)+0.405837E1.*cos(12.*t)+(-0.4004E1).*cos(13.*t)+( ...
%   -0.187109E1).*cos(14.*t)+(-0.513127E0).*cos(15.*t)+(-0.223154E1).* ...
%   cos(16.*t)+(-0.299107E1).*cos(17.*t)+(-0.179691E1).*cos(18.*t)+ ...
%   0.677947E0.*cos(19.*t)+(-0.778444E0).*cos(20.*t)+(-0.194499E3).* ...
%   sin(t)+0.548228E2.*sin(2.*t)+(-0.145397E2).*sin(3.*t)+0.336624E0.* ...
%   sin(4.*t)+(-0.253067E1).*sin(5.*t)+0.262837E0.*sin(6.*t)+ ...
%   0.388591E1.*sin(7.*t)+0.437869E1.*sin(8.*t)+(-0.201029E1).*sin(9.* ...
%   t)+(-0.363847E1).*sin(10.*t)+(-0.342847E1).*sin(11.*t)+( ...
%   -0.269287E1).*sin(12.*t)+(-0.34024E1).*sin(13.*t)+(-0.330055E1).* ...
%   sin(14.*t)+0.101505E1.*sin(15.*t)+0.954781E0.*sin(16.*t)+( ...
%   -0.499421E0).*sin(17.*t)+0.178974E1.*sin(18.*t)+0.166851E1.*sin( ...
%   19.*t)+(-0.110728E1).*sin(20.*t);
% Curve_y=0.146617E3+(-0.100911E3).*cos(t)+0.43755E2.*cos(2.*t)+( ...
%   -0.335871E2).*cos(3.*t)+(-0.723942E1).*cos(4.*t)+0.100098E2.*cos( ...
%   5.*t)+0.347639E1.*cos(6.*t)+(-0.673488E1).*cos(7.*t)+(-0.260574E2) ...
%   .*cos(8.*t)+(-0.113183E2).*cos(9.*t)+0.420296E1.*cos(10.*t)+( ...
%   -0.336902E1).*cos(11.*t)+0.646439E0.*cos(12.*t)+(-0.456828E1).* ...
%   cos(13.*t)+(-0.487865E1).*cos(14.*t)+0.312336E1.*cos(15.*t)+ ...
%   0.182675E1.*cos(16.*t)+(-0.101101E1).*cos(17.*t)+0.849381E0.*cos( ...
%   18.*t)+0.107306E1.*cos(19.*t)+(-0.111791E1).*cos(20.*t)+ ...
%   0.394822E2.*sin(t)+(-0.180802E2).*sin(2.*t)+0.716137E1.*sin(3.*t)+ ...
%   (-0.386693E2).*sin(4.*t)+(-0.849595E1).*sin(5.*t)+0.273687E2.*sin( ...
%   6.*t)+(-0.368229E1).*sin(7.*t)+0.846716E1.*sin(8.*t)+0.154497E2.* ...
%   sin(9.*t)+0.171189E1.*sin(10.*t)+(-0.909299E0).*sin(11.*t)+( ...
%   -0.250452E0).*sin(12.*t)+(-0.338652E1).*sin(13.*t)+(-0.324023E1).* ...
%   sin(14.*t)+0.108961E1.*sin(15.*t)+0.276616E-1.*sin(16.*t)+( ...
%   -0.149295E1).*sin(17.*t)+(-0.43546E0).*sin(18.*t)+(-0.989119E0).* ...
%   sin(19.*t)+0.925713E0.*sin(20.*t);
% Curve_x=R0.*(Curve_x/max(abs(Curve_x(:))));
% Curve_y=R0.*(Curve_y/max(abs(Curve_y(:))));

% 3. Expression of GPOV beams
Y = zeros ( N , M ) ;
dCurve_x = gradient ( Curve_x , T / ( P - 1 ) ) ;
dCurve_y = gradient ( Curve_y , T / ( P - 1 ) ) ;
Jacobiant = Curve_x .* dCurve_y - Curve_y .* dCurve_x ;
% Two integral is not permitted in matlab 
% (But cumtrapz works due to the pure numerical method)
value1 = cumtrapz ( t , Jacobiant ) ;
l = 3 ; % Topological Charge
sigma = l * 2 * pi / trapz ( t , Jacobiant ) ;
lamda = 632.8e-9 ; % wavelength
f = 0.2 ; % focal length
% GPOV
for i = 1 : N
    for j = 1 : M
%         xx = x ( i , j , 1 ) * ( M - 1 ) / 4 / a^2 * lamda * f ; % Frequency Domain
%         yy = y ( i , j , 1 ) * ( N - 1 ) / 4 / b^2 * lamda * f ; % Frequency Domain
        xx = x ( i , j , 1 ) ;
        yy = y ( i , j , 1 ) ;
        Phi = exp ( 1i * 2 * pi / lamda / f * ( xx .* Curve_x + yy .* Curve_y )...
            + 1i * sigma * value1 ) ;
        % Y
        value2 = Phi .* Jacobiant ;
        Y ( i , j ) = trapz ( t , value2 ) / lamda / f ;    
    end
end
Y = Y / max ( max ( abs ( Y ) ) ) ;
% C
C = fftshift ( fft2 ( Y ) ) ; % May induce a shift bug in BPM (so we use FGPOV for masks).
C = C / max ( max ( abs ( C ) ) ) ;
end