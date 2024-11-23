function [fval] = tima_fval_cvm(observed,modeled,~,~)
%TIMA_FVAL_CVM
%   function to calculate Kolmogorov-Smirnov Test Statistic
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

   % Sort the observations and modeled values
    sortedObs = sort(observed);
    sortedMod = sort(modeled);
    
    % Compute the empirical CDF of the observed data
    n = length(observed);
    F_n = (1:n) / n; % Empirical CDF values
    x_empirical = sortedObs; % Corresponding sorted observed values
    
    % Compute the theoretical CDF of the modeled data
    F_m = (1:n) / n; % CDF of modeled values (assuming uniformly distributed)
    x_model = sortedMod; % Corresponding sorted modeled values

    % Compute the areas under the curves
    F_n_interp = interp1(x_empirical, F_n, sortedObs, 'linear', 'extrap');
    F_m_interp = interp1(x_model, F_m, sortedObs, 'linear', 'extrap');
    
    % Calculate the Cramér-von Mises statistic
    fval = sum((F_n_interp - F_m_interp).^2) / n;
end

