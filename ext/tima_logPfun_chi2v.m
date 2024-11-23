function [logfval] = tima_logPfun_chi2v(observed,modeled,error)
%TIMA_LOGPFUN_CHI2V
%   function to calculate reduced chi squared goodness of fit
%
% Usage in TIMA
%   logfval = @(theta) tima_logPfun_chi2v(TData.temps_to_fit(MData.fit_ind),tima_formod_subset(theta,MData.fit_ind,formod),MData.erf(TData.temps_to_fit(MData.fit_ind)))
%   logfval = @(theta) tima_logPfun_chi2v([TData.temps_to_fit(MData.fit_ind);TData.temps_to_fit_II(MData.fit_ind)],[tima_formod_subset(theta,MData.fit_ind,formod);tima_formod_subset(theta,MData.fit_ind,formod_II)],[MData.erf(TData.temps_to_fit(MData.fit_ind); MData.erf(TData.temps_to_fit_II(MData.fit_ind))]]))
%
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

logfval = -0.5*sum((observed-modeled).^2./(error).^2);% a cell of function handles returning the log probality of a each outcom
end

