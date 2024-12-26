function [fval] = tima_fval_DTW(observed,modeled,error,nvars)
%TIMA_FVAL_DTW
%   function to calculate modified reduced chi squared goodness of fit
%   accounting for distance time warping (DTW)
%   Chi^2_v reduced chi squared = goodness of fit
%       -chi^2_v  1 Fit is poor. The chosen model is unlikely to fit the data.
%       -chi^2_v  1 Fit does not fully capture the data. It's possible there are too few data points or the uncertainties are underestimated.
%       -chi^2_v  1 Fit is too good (overfitting). The model fits noise in the data or the uncertainties are overestimated (can be caused by having too many free parameters in fitting).
%       In general chi^2_v ~= 1 means that the fit is doing a good job without overfitting. While these guidelines don't replace full statistical tests, they can be quick check of goodness of fit.
%
% Usage in TIMA
%   Obj = @(theta) tima_fval_DTW(Temps_Obs,tima_formod_subset(theta,UAV_flight_ind,formod),MData.erf(Temps_Obs),MData.nvars)

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
[~,io,im] = dtw(observed,modeled);
fval = sum((observed(io)-modeled(im)).^2./error(io).^2)/(length(io)-nvars);
end

