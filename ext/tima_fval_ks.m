function [fval] = tima_fval_ks(observed,modeled,~,~)
%TIMA_FVAL_KS
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

    % Set number of bootstrap samples
    num_samples = 1000;
    W2_samples = zeros(num_samples, 1);
    
    % Get length of observed data
    n = length(observed);

    for i = 1:num_samples
        % Generate resampled observed data based on uncertainty
        resampled_obs = observed + uncertainty .* randn(size(observed));
        
        % Sort the resampled observations and modeled values
        sortedObs = sort(resampled_obs);
        sortedMod = sort(modeled);
        
        % Compute the empirical CDF of the resampled data
        F_n = (1:n) / n; % Empirical CDF values
        x_empirical = sortedObs; % Corresponding sorted observed values
        
        % Compute the theoretical CDF of the modeled data
        F_m = (1:n) / n; % CDF values for modeled data (assumed uniform again)
        x_model = sortedMod; % Corresponding sorted modeled values

        % Interpolating the empirical CDF and modeled CDF
        F_n_interp = interp1(x_empirical, F_n, sortedObs, 'linear', 'extrap');
        F_m_interp = interp1(x_model, F_m, sortedObs, 'linear', 'extrap');

        % Calculate the Cramér-von Mises statistic
        W2_samples(i) = sum((F_n_interp - F_m_interp).^2) / n;
    end
    
    % Calculate the final Cramér-von Mises statistic for the original observed data
    % Using average over bootstrap samples. This may not be strictly necessary,
    % but provides information on the variability of W2.
    W2 = mean(W2_samples);
    
    % Display the final result and variability
    fprintf('Cramér-von Mises statistic (considering uncertainty): W2 = %.6f\n', W2);
end

