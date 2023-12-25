    function [mu, yi] = jac2Plt (m, ab, damping, xi) 
%%  function [mu, yi] = jac2Plt (m, ab, damping, xi) 
%% computes pn and pn (xi) where pn is a Cheb. epansion of
%% a step function with value = 1 in [a b] and zero
%% outside [a b].
% By Shashanka Ubaru and Yousef Saad
%% m       = degree
%% ab      = ab(1) = a, ab(2) = b  with -1 <= a <= b <= 1  
%%    a = -1 --> high-pass
%%    b =  1 -->low-pass
%%    else   --> midpass [barrier]
%% damping = real number indicating whether or not to add damping
%%    damping == 1 --> Jackson, Damping == 0 --> Standard Cheb, exp
%% xi      = (optional) a vector of values where p(xi) is to be computed.
%% Note xi's must be in [-1 1] x
%%      
%% return : 
%% mu      =  expansion coefficients (in cheb. pols) coefs.
%% yi      =  pn(xi)   yi =0 if no points xi provided
%%
 alpha1 = ab(1);
 alpha2 = ab(2);
%%-------------------- expansion coefficients
 thetJ = pi/(m+2); 
 thetL = pi/(m+1); 
 a1 = 1/(m+2);
 a2 = sin(thetJ);
 beta1 = acos(alpha1);
 beta2 = acos(alpha2);
for k=0:m 
%%-------------------- generate Jackson coefficients
 if (damping == 0) 
     jac = 1;
 elseif (damping  == 1) 
     %%-------------------- Note: slightly simpler formulas for
     %%jackson: 
     jac = a1*sin((k+1)*thetJ)/a2 + (1-(k+1)*a1)*cos(k*thetJ) ; 
     %%-------------------- Lanczos sigma-damping: 
 elseif (damping == 2)  
     jac = 1;
     if (k>0)
         jac = sin(k*thetL)/(k*thetL);
     end
 end
 if (k == 0)  
     mu(k+1) = -jac*(beta2-beta1)/pi;   
 else
     mu(k+1) = -2*jac*(sin(k*beta2)-sin(k*beta1))/(pi*k);  
 end 
end
%%-------------------- return if no plotting to be done
   if (nargin < 4)
       yi = 0;
       return 
   end
   
%%-------------------- compute p(xi)
   n = size(xi,1) ; 
   yi = zeros(n,1);
   vkm1 = zeros(n,1);
   vk   = ones(n,1);
   for k=0:m
%%-------------------- accumulation of vector.
       yi = yi + mu(k+1)*vk;
       scal = 2;
       if (k == 0) , scal = 1;, end;
       vkp1 = scal .*xi .* vk - vkm1;
       vkm1 = vk; 
       vk = vkp1;
 end
