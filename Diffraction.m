%%  Beam Propagation Method
function [ E , Ezi , Exoz ] = Diffraction ( E0 , BPM_Method )
% 1. Environrment setup
% Global invariants
global a;global b;global c;global M;global N;global MN;
global xa;global ya;global za;global x;global y;global z;global r;
lamda = 632.8e-9 ;
% Method Options (0 for D-FFT & 1 for S-FFT) 
if nargin<2
    BPM_Method = 0;
end
if isempty(E0)==1
    error('No beam input'); 
end
if isempty(BPM_Method)==1
   BPM_Method=0;
end
% 2. Diffraction 
E=zeros(N,M,MN);
% G0=fftshift(fft2(ifftshift(E0)));
G0=fftshift(fft2(E0));
% Initial plane excepted
k0=1+floor(MN/2);
E(:,:,k0)=E0;
% Beam propagation method (designed for odd number)
% Negetive Position calculation
for k=1:floor(MN/2)
        if BPM_Method == 0
    % D-FFT method (caution for opposite calculation due to the huge value)
    Gz_conj=conj(G0).*exp(1i*2*pi*abs(z(1,1,k)).*sqrt(1/(lamda.^2)-...
        (x(:,:,1)/a/a/4*(M-1)).^2-(y(:,:,1)/b/b/4*(N-1)).^2));
%     Gz=G0.*exp(1i*2*pi*z(1,1,k)/lamda.*(1-lamda^2/2*((x/4/a/a*M).^2+(y/4/b/b*N).^2)));
%     Ez=fftshift(ifft2(ifftshift(Gz)));
    % Circ Filter
    Gz=conj(Gz_conj);
    % When using indicators, circle is fine; Using filter, ellipse.
    % Not necessary due to the evanescent wave 
    T=CircFilter(1/lamda*4*a*a/(M-1),1/lamda*4*b*b/(N-1)); 
    Gz=Gz.*T;
    Ez=ifft2(ifftshift(Gz));
        else
    % S-FFT method
    a2=(M-1)*lamda*abs(z(1,1,k))/4/a;b2=(N-1)*lamda*abs(z(1,1,k))/4/b;
    xa2=linspace(-a2,a2,M);ya2=linspace(b2,-b2,N);
    [x2,y2]=meshgrid(xa2,ya2);r2=sqrt(x2.^2+y2.^2);
    Ez2=fftshift(fft2(E0.*exp(1i*pi/lamda/z(1,1,k)*r(:,:,1).^2)));
    % Change coordinate
    Ez=exp(1i*2*pi/lamda*z(1,1,k))/1i/lamda/z(1,1,k)...
        *exp(1i*pi/lamda/z(1,1,k)*r2.^2).*Ez2;
        end   
%     % Normalized ( Generally unnecessary)
%     Ez=Ez./max(abs(Ez(:)));
    for i=1:N
        for j=1:M
% Fill Ez into Ezz
          E(i,j,k)=Ez(i,j);     
        end
    end
end
% Positive calculation
for k=(floor(MN/2)+2):MN
        if BPM_Method == 0
    % D-FFT method (caution for opposite calculation due to the huge value)
    Gz=G0.*exp(1i*2*pi*abs(z(1,1,k)).*sqrt(1/(lamda.^2)-...
        (x(:,:,1)/a/a/4*(M-1)).^2-(y(:,:,1)/b/b/4*(N-1)).^2));
%     Gz=G0.*exp(1i*2*pi*z(1,1,k)/lamda.*(1-lamda^2/2*((x/4/a/a*M).^2+(y/4/b/b*N).^2)));
%     Ez=fftshift(ifft2(ifftshift(Gz)));
    % Circ Filter
    T=CircFilter(1/lamda*4*a*a/(M-1),1/lamda*4*b*b/(N-1));
    Gz=Gz.*T;
    Ez=ifft2(ifftshift(Gz));
        else
    % S-FFT method
    a2=(M-1)*lamda*abs(z(1,1,k))/4/a;b2=(N-1)*lamda*abs(z(1,1,k))/4/b;
    xa2=linspace(-a2,a2,M);ya2=linspace(b2,-b2,N);
    [x2,y2]=meshgrid(xa2,ya2);r2=sqrt(x2.^2+y2.^2);
    Ez2=fftshift(fft2(E0.*exp(1i*pi/lamda/z(1,1,k)*r(:,:,1).^2)));
    % Change coordinate
    Ez=exp(1i*2*pi/lamda*z(1,1,k))/1i/lamda/z(1,1,k)...
        *exp(1i*pi/lamda/z(1,1,k)*r2.^2).*Ez2;
        end   
%     % Normalized
%     Ez=Ez./max(abs(Ez(:)));
    for i=1:N
        for j=1:M
% Fill Ez into Ezz
          E(i,j,k)=Ez(i,j);     
        end
    end
end
% 3.Ezi
Ezi=E(:,:,MN);
% 4.xoz 
Dxoz=floor(N/2);
Exoz=zeros(M,MN);
for j=1:M
    for k=1:MN
        Exoz(j,k)=E(Dxoz,j,k);
    end
end
end