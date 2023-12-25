function yi = TcPlt (mu, xi) 
%% function yi = TcPlt (mu, xi) 
%% 
%% Plot the Chebyshev polynomial.  Jackson damping is not used here but
%% is assumed to be multiplied by mu outside this routine.
n = size(xi,1); 
m = length(mu);
%%-------------------- compute p(xi)
yi = zeros(n,1);
vkm1 = zeros(n,1);
vk   = ones(n,1);
for k=0:m-1
	%%-------------------- accumulation of vector.
	yi = yi + mu(k+1)*vk;
	scal = 2;
	if (k == 0) 
		scal = 1;
	end
	vkp1 = scal .*xi .* vk - vkm1;
	vkm1 = vk; 
	vk = vkp1;
end
yi = yi ./ sqrt(1 - xi .^2);
