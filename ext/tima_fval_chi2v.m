function [fval] = tima_fval_chi2v(observed,modeled,error,nvars)
%TIMA_FVAL_CHI2V
%   function to calculate reduced chi squared goodness of fit
%
% Usage in TIMA
%   Obj = @(theta) tima_fval_chi2v(Temps_Obs,tima_formod_subset(theta,UAV_flight_ind,formod),MData.erf(Temps_Obs),MData.nvars)

% Inputs
%   observed: array of observed values for fitting (vector)
%   modeled: array of modeled vlaues the same size as observed (vector)
%   error: array of uncertainty values the same size as observed (vector)
%   nvars: number of free variable (vector)
%
% Outputs
%   Output: fit value (1=good fit,<1=overfit,>1=imperfect fit)
%
% Author
%    Ari Koeppel, 2021

fval = sum((observed-modeled).^2./error.^2)/(length(observed)-nvars);
end

