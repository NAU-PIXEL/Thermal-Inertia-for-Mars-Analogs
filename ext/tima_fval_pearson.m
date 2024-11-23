function [fval] = tima_logPfun_pearson(inputArg1,inputArg2)
%TIMA_LOGPFUN_PEARSON Summary of this function goes here
%   Detailed explanation goes here
    % Sort the observed and modeled data
    observed_sorted = sort(observed);
    modeled_sorted = sort(modeled);

    % Calculate the CDF for the observed and modeled arrays
    n = numel(observed);  % Number of observed data points
    m = numel(modeled);   % Number of modeled data points
    cdf_observed = (1:n)' / n;  % Empirical CDF of observed
    cdf_modeled = (1:m)' / m;    % Empirical CDF of modeled

    % Combine values for CDF calculation
    all_values = unique([observed_sorted; modeled_sorted]);
    cdf_obs = zeros(size(all_values));
    cdf_mod = zeros(size(all_values));

    % Calculate the empirical CDF for the observed data
    for i = 1:length(all_values)
        cdf_obs(i) = sum(observed_sorted <= all_values(i)) / n;
    end

    % Calculate the empirical CDF for the modeled data
    for i = 1:length(all_values)
        cdf_mod(i) = sum(modeled_sorted <= all_values(i)) / m;
    end

    % Calculate the KS statistic
    ks_statistic = max(abs(cdf_obs - cdf_mod));
end

