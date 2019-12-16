%Computes the lognormal fit parameters given an input distribution
%Gabriel Martine
%Last updated: May 13 2019

function [mu, sigma, err_mu, err_sigma] = LognormalParams(dist)
	[parmhat, parmcil] = lognfit(dist);
	mu = parmhat(1); err_mu = abs(parmcil(2) - mu); 
	sigma = parmhat(2); err_sigma = abs(parmcil(4) - sigma);
end
