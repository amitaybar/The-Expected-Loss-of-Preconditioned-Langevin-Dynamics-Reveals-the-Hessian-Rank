function [r,eps,zz,z1,xx,yy]=KPM_DOS_and_Rank(A,lmin,lmax,deg,nvecs,Npts,AutoGap,n)
%% function [r,g,zz,z1,xx,yy]=KPM_DOS_and_Rank(A,lmin,lmax,deg,nvecs,Npts,AutoGap)
% computes the DOS, estimates tthe threshold \eps 
% and then computes the eigencount above this threshold (the rank)
% using KPM and Chebyshev polynomial filters
%By Shashanka Ubaru and Yousef Saad

%Inputs:
% A= input matrix, [lmin,lmax] =min and max eigenvalues of A
% deg = degree of Chebyshev polynomial
% nvecs = number of samples
% Npts = number of points in DOS
% AutoGap = flag to estimate DOS and the gap

% Outputs:
% r = estimated rank
% g = threshold 
% z1 = rank for different starting vectors
% zz = cumulative average of ranks
% xx = x axis of DOS plot
% yy = Values of DOS plot 

%% intialization
% n = size(A,2);
dd = (lmax - lmin)/2;
cc = (lmax + lmin)/2;
B = @(v)(A(v) - cc*v)/dd;
z1=zeros(nvecs,1);
zz=zeros(nvecs,1);

%% Form v'Tk(A)v for nvecs starting vectors
Y=zeros(deg+1,nvecs);
for ii = 1:nvecs
    v = randn(n,1) ;
    v = v/norm(v,2);
    z= Chebyshev_vTAv(B, v, deg);
    Y(:,ii)=z;
end

%% Find DOS and Threshold
if(AutoGap)
    damping  = 1; % type of damping: 0 - no damping; 1 - Jackson and 2 - Lanczos
    [xx, yy] = kpmchebdos_new(Y,Npts,lmin,lmax,damping);
    %%-------------------- locate gap
    eps = find_gap(xx,yy,n);
else
    xx=zeros(Npts,1);
    yy=zeros(Npts,1);
    eps=0.5;  %set a threshold 
end

%% Find approximate rank
ab0 = [eps , lmax];
ab = (ab0 - [cc cc] ) ./ dd;
cnt_est = 0;
[mu,~] = jac2Plt(deg, ab, 2) ;

for l = 1:nvecs
    vk=Y(:,l);
    t=sum(mu'.*vk);
    z1(l) = t*n;
    cnt_est = (cnt_est+t);
    zz(l) = n*(cnt_est/l);
end
r = zz(nvecs);
end

