%% CircFilter
function T=CircFilter(a1,a2)
%%  0.
% 1.Global invariants
global a;global b;global c;global M;global N;global MN;
global xa;global ya;global za;global x;global y;global z;global r;
% 2. Parameter
if nargin<2
    a1=a;a2=b;
end
if isempty(a1)==1
    a1=a;
end
if isempty(a2)==1
    a2=b;
end

%%  1.
T=zeros(N,M);
for i=1:N
    for j=1:M
        if (((j-(M-1)/2-1)*2*a/(M-1))^2/(a1^2)+((i-(N-1)/2-1)*2*b/(N-1))^2/(a2^2))<1
            T(i,j)=1;
        end
    end
end
end