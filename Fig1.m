%  Figure1
%
% %%%%%%%%%%%%%%%%%%%%%         Description        %%%%%%%%%%%%%%%%%%
% 1. Visualize the intensity and phase distribution of simulated beams.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by: Yue chen
% Department of Modern Physics, University of Science and Technology of
% China (USTC)
% December. 4, 2018
% More imformation£ºchenyue8@mail.ustc.edu.cn

function Fig1 ( E10 , Ezi1 , Exoz1 , E20 , Ezi2 , Exoz2 )
global a;global b;global c;global M;global N;global MN;
global xa;global ya;global za;global x;global y;global z;global r;
lamda = 632.8e-9 ;
f = 0.2 ;

% 
figure ( 1 )
subplot ( 2 , 5 , 1 )
imagesc ( abs ( E10 ) .^ 2 ) ; colormap ( gray ) ;
title ( 'C z=0 Intensity' ) ;
xlabel ( 'x/m' ) ;
ylabel ( 'y/m' ) ;
set ( gca , 'xtick' , linspace ( 1 , M , 5 ) , 'xticklabel' , sprintfc (...
    '%g' , linspace ( -a , a , 5 ) ) ,...
    'ytick' , linspace ( 1 , N , 5 ) , 'yticklabel' , sprintfc( '%g' , ...
    linspace ( b , -b , 5 ) ) ) ;
subplot ( 2 , 5 , 2 )
imagesc ( mod ( angle ( E10 ) , 2 * pi ) ) ; colormap ( gray ) ;
title ( 'C z=0 Phase' ) ;
xlabel ( 'x/m' ) ;
ylabel ( 'y/m' ) ;
set ( gca , 'xtick' , linspace ( 1 , M , 5 ) , 'xticklabel' , sprintfc (...
    '%g' , linspace ( -a , a , 5 ) ) ,...
    'ytick' , linspace ( 1 , N , 5 ) , 'yticklabel' , sprintfc ( '%g' ,...
    linspace ( b , -b , 5 ) ) ) ;
subplot ( 2 , 5 , 3 )
imagesc ( abs ( Ezi1 ) .^ 2 ) ; colormap ( gray ) ;
title ( ['C Intensity at z=',num2str(c)] ) ;
xlabel ( 'x/m' ) ;
ylabel ( 'y/m' ) ;
set ( gca , 'xtick' , linspace ( 1 , M , 5 ) , 'xticklabel' , sprintfc (...
    '%g' , linspace ( -a , a , 5 ) ) ,...
    'ytick' , linspace ( 1 , N , 5 ) , 'yticklabel' , sprintfc ( '%g' , ...
    linspace ( b , -b , 5 ) ) ) ;
subplot ( 2 , 5 , 4 )
imagesc ( mod ( angle ( Ezi1 ) , 2 * pi ) ) ; colormap ( gray ) ;
title ( ['C Phase at z=',num2str(c)] ) ;
xlabel ( 'x/m' ) ;
ylabel ( 'y/m' ) ;
set ( gca , 'xtick' , linspace ( 1 , M , 5 ) , 'xticklabel' , sprintfc(...
    '%g' , linspace ( -a , a , 5 ) ) ,...
    'ytick' , linspace ( 1 , N , 5 ) , 'yticklabel' , sprintfc ( '%g' ,...
    linspace(b,-b,5)));
subplot ( 2 , 5 , 5 )
imagesc ( abs ( Exoz1 ) .^ 2 ) ; colormap ( gray ) ;
title ( 'C y=0 Intensity' ) ;
xlabel ( 'z/m' ) ;
ylabel ( 'x/m' ) ;
set ( gca , 'xtick' , linspace ( 1 , MN , 5 ) , 'xticklabel' , sprintfc(...
    '%g' , linspace ( -c , c , 5 ) ) ,...
    'ytick' , linspace ( 1 , M , 5 ) , 'yticklabel' , sprintfc ( '%g' , ...
    linspace ( a , -a , 5 ) ) ) ;
subplot(2,5,6)
imagesc(abs(E20).^2);colormap(gray);
title('Y z=0 Intensity');
xlabel('x/m');
ylabel('y/m');
set(gca,'xtick',linspace(1,M,5),'xticklabel',sprintfc('%g',linspace(-(M-1)/4/a* lamda * f,(M-1)/4/a* lamda * f,5)),...
'ytick',linspace(1,N,5),'yticklabel',sprintfc('%g',linspace((N-1)/4/b* lamda * f,-(N-1)/4/b* lamda * f,5)));
subplot(2,5,7)
imagesc(mod(angle(E20),2*pi));colormap(gray);
title('Y z=0 Phase');
xlabel('x/m');
ylabel('y/m');
set(gca,'xtick',linspace(1,M,5),'xticklabel',sprintfc('%g',linspace(-(M-1)/4/a* lamda * f,(M-1)/4/a* lamda * f,5)),...
'ytick',linspace(1,N,5),'yticklabel',sprintfc('%g',linspace((N-1)/4/b* lamda * f,-(N-1)/4/b* lamda * f,5)));
subplot(2,5,8)
imagesc(abs(Ezi2).^2);colormap(gray);
title(['Y Intensity at z=',num2str(c)]);
xlabel('x/m');
ylabel('y/m');
set(gca,'xtick',linspace(1,M,5),'xticklabel',sprintfc('%g',linspace(-(M-1)/4/a* lamda * f,(M-1)/4/a* lamda * f,5)),...
'ytick',linspace(1,N,5),'yticklabel',sprintfc('%g',linspace((N-1)/4/b* lamda * f,-(N-1)/4/b* lamda * f,5)));
subplot(2,5,9)
imagesc(mod(angle(Ezi2),2*pi));colormap(gray);
title(['Y Phase at z=',num2str(c)]);
xlabel('x/m');
ylabel('y/m');
set(gca,'xtick',linspace(1,M,5),'xticklabel',sprintfc('%g',linspace(-(M-1)/4/a* lamda * f,(M-1)/4/a* lamda * f,5)),...
'ytick',linspace(1,N,5),'yticklabel',sprintfc('%g',linspace((N-1)/4/b* lamda * f,-(N-1)/4/b* lamda * f,5)));
subplot(2,5,10)
imagesc(abs(Exoz2).^2);colormap(gray);
title('Y y=0 Intensity');
xlabel('z/m');
ylabel('x/m');
set(gca,'xtick',linspace(1,MN,5),'xticklabel',sprintfc('%g',linspace(-c,c,5)),...
'ytick',linspace(1,M,5),'yticklabel',sprintfc('%g',linspace(-(M-1)/4/a* lamda * f,(M-1)/4/a* lamda * f,5)));
% figure(6+4*yanshe)
% imagesc(mod(Pxoz,2*pi));colormap(gray);
end