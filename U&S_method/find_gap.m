function [gap] = find_gap(xx,yy,n)
%% function [gap] = find_gap(xx,yy,n)
%% finds first gap in DOS curves..
%By Shashanka Ubaru and Yousef Saad
%%
npts = length(yy);
%%-------------------- tolerance for high values
tol1 = (npts+1)*0.5;
%%-------------------- smooth the curve
p = 3;   % p = 0 ----> no smoothing
zz = abs(yy);
for ii=2:p
    zz(1:npts-ii) = zz(1:npts-ii) + yy(ii+1:npts);
end
yy = zz/(p+1) ;
%%-------------------- derivatives are often better
%%                     at spotting gap:
tt = yy(2:npts)-yy(1:npts-1);  % * (xx(2)-xx(1));
tt(npts) = 0;
%%-------------------- get 1st point of increase
ysum = 0;
for i=1:npts
    %%-------------------- spot first time derivative turns>0
    ysum = ysum+yy(i);
    %%-------------------- discard high values
    if (yy(i)>tol1)
        continue
    end
    if (tt(i)>=-0.01)
        indx1 = i;
        break;
    end
end
ysum = ysum/(npts+1);
disp('ysum indx1')
disp([ysum, indx1])
%%-------------------- spot first time yy becomes significant again
tol  = (1 - ysum)*0.25
indx2=[];
for i=indx1+1:npts
    %%-------------------- spot first time derivative becomes>0
    if (yy(i) >= tol)
        indx2 = i;
        break;
    end
end
if(isempty(indx2))
    indx2=indx1+1;
end
disp('indx1 indx2 = ')
[indx1  indx2]
%%
[yy(indx1), yy(indx2), xx(indx1), xx(indx2)]
%%-------------------- middle point between indx and indx2
%%%gap = xx(indx1);
gap = (xx(indx1)+xx(indx2))/2;
