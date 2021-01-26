%%  Figure2
function Fig2(DMDpixels1,DMDpixels2,WorN)
global a;global b;global c;global M;global N;global MN;
global xa;global ya;global za;global x;global y;global z;global r;
figure(2)
if isempty(WorN)==1
    WorN=1;
end
% % fill up to 1024*768
% DMDpixels11=zeros(768,1024);
% DMDpixels22=zeros(768,1024);
% DMDpixelsG11=zeros(768,1024);
% DMDpixelsG22=zeros(768,1024);
% DMDpixels11(:,128:895)=DMDpixels1;
% DMDpixels22(:,128:895)=DMDpixels2;
% DMDpixelsG11(:,128:895)=DMDpixelsG1;
% DMDpixelsG22(:,128:895)=DMDpixelsG2;

% % 768*768 fill up to 1920*1080 (gray)
% DMDpixels11=zeros(1080,1920);
% DMDpixels22=zeros(1080,1920);
% DMDpixels11(157:924,577:1344)=abs(DMDpixels1);
% DMDpixels22(157:924,577:1344)=abs(DMDpixels2);

% 1080*1080 fill up to 1920*1080 (gray)
DMDpixels11=zeros(1080,1920);
DMDpixels22=zeros(1080,1920);
DMDpixels11(:,421:1500)=abs(DMDpixels1);
DMDpixels22(:,421:1500)=abs(DMDpixels2);

% % 768*768 fill up to 1920*1080 (logical)
% DMDpixels11=false(1080,1920);
% DMDpixels22=false(1080,1920);
% DMDpixels11(157:924,577:1344)=DMDpixels1;
% DMDpixels22(157:924,577:1344)=DMDpixels2;

% % 1080*1080 fill up to 1920*1080 (logical)
% DMDpixels11=false(1080,1920);
% DMDpixels22=false(1080,1920);
% DMDpixels11(:,421:1500)=DMDpixels1;
% DMDpixels22(:,421:1500)=DMDpixels2;

% image and write
subplot(1,2,1)
imagesc(DMDpixels11);colormap(gray);
title('C Hologram');
subplot(1,2,2)
imagesc(DMDpixels22);colormap(gray);
title('Y Hologram');
if WorN==1
    imwrite(DMDpixels11,'.\[6]HologramGeneration\Hologram_C.bmp');
    imwrite(DMDpixels22,'.\[6]HologramGeneration\Hologram_Y.bmp');
end
end
