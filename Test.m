figure;
imagesc(abs(Y0).^2);colormap('hot');axis off;axis square;
figure;
imagesc(mod(angle(Y0),2*pi));colormap('gray');axis off;axis square;