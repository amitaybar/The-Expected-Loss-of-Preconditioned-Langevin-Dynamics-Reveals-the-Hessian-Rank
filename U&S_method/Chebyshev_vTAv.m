    function y = Chebyshev_vTAv(A, v, deg)
%%    function y = Chebyshev_vTAv(A, v, deg)
%% 
%% computes  y(k) = v'Tk(A)v, where Tk is a Cheb. polynomial 

 n = size(v,1) ; 
 v = v(:);
%%-------------------- compute pn(A)v
   y = zeros(deg,1);
   vkm1 = zeros(n,1);
%%         T0(xi) = ones*v = v
   vk   = v;
   y(1,1)=1;
   for k=1:deg
%%-------------------- accumulation of vector.
        scal = 2;
        if (k == 1) 
            scal = 1;
        end
        vkp1 = scal* (A(vk)) - vkm1;
        vkm1 = vk; 
        vk = vkp1;
        y(k+1,1) = (v'*vk);
   end
 end
