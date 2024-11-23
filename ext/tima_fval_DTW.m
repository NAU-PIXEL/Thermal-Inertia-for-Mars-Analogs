function [fval] = tima_fval_DTW(observed,modeled,~,~)
%TIMA_FVAL_DTW
%   function to calculate goodness of fit by distance time warp
%
% Inputs
%   observed: array of observed values for fitting (vector)
%   modeled: array of modeled vlaues the same size as observed (vector)
%   error: array of uncertainty values the same size as observed (vector)
%   nvars: number of free variable (vector)
%
% Outputs
%   Output: fit value
%
% Author
%    Ari Koeppel, 2021

[fval] = dtw(observed, modeled);
end

