function [xx, yy] = kpmchebdos_new(Y, Npts, lm, lM, damping)
%% function [xx, yy] = kpmchebdos_new(Y, Npts, lm, lM, damping)
% Computes the DOS of the matrix using KPM
% By Shashanka Ubaru and Yousef Saad

%Inputs:
% Y=  matrix v'Tk(A)v for m degree and nvecs starting vectors
%[lm,lM] =min and max eigenvalues of A
% Npts = number of points in DOS
%  damping =  type of damping: 
% 0 - no damping; 1 - Jackson and 2 - Lanczos
   

% Outputs:
% xx = x axis of DOS plot
% yy = Values of DOS plot 

%% --------------------
h = 2/(Npts+1);
pts = [-1+h:h*0.5:0,0+h:2*h:1-h]';
%pts = [-1+h:h:1-h]';
[deg,nvec]=size(Y);

ctr = (lm+lM)/2;
wid = (lM-lm)/2 + 2*eps;

mu = zeros(deg,1);
for m = 1 : nvec
    wk   = Y(:,m);
    thetJ = pi/(deg+2);
    thetL = pi/(deg+1);
    a1 = 1/(deg+2);
    a2 = sin(thetJ);
    for k=0:deg-1
        if (damping == 0)
            jac = 1;
        elseif (damping  == 1)
            %%----- Note: slightly simpler formulas for jackson:
            jac = a1*sin((k+1)*thetJ)/a2+(1-(k+1)*a1)*cos(k*thetJ);
            %%-------------------- Lanczos sigma-damping:
        elseif (damping == 2)
            jac = 1;
            if (k>0)
                jac = sin(k*thetL)/(k*thetL);
            end
        end
        mu(k+1) = mu(k+1)+jac*(wk(k+1));
    end
end

%%-------------------- Plot
mu = 2*mu ./ (nvec*pi);  %% scaling in formula
mu(1) = mu(1)/2;         %% first term is different

xx  = ctr + wid*pts;
y2 = TcPlt(mu, pts);
yy  = y2';
%-------------------- yy = (nptsC/nptsR) *yy /sum(yy);
yy  = yy / (sum(yy)*(xx(2)-xx(1)));
