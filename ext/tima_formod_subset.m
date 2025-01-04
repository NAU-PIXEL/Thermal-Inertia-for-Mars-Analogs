function [Output] = tima_formod_subset(FitParams,indecies,formod)
%% TIMA_FORMOD_SUBSET
%   function to run forwad model but only retrieve results for a specific subset of indecies
%
% Syntax
%   [Output] = tima_formod_subset(FitParams,indecies,formod)

% Inputs
%   FitParams: parameters for fitting (vector)
%   indecies: indecies of full dataset for which observational data exists
%       (vector)
%   formod: forward model function (function handle)
%
% Outputs
%   Output: forward model output constrained to specific indecies
%
% Author
%    Ari Koeppel, 2021

Results = formod(FitParams);
Output = Results(indecies);
end
