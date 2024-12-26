function [fval] = tima_fval_dchi2v_dt(observed,modeled,error,nvars)
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

    dobs = gradient(observed);
    dmod = gradient(modeled);

    % Remove the last element to match sizes (ddata1 and ddata2 will be 1 shorter)
    dobs = movmean(dobs,60);
    dmod = movmean(dmod,60);
    unc = error;

    % Compute the chi-squared statistic
    chi_squared = sum(((dobs - dmod).^2) ./ unc.^2);

    % Degrees of freedom
    nu = length(dobs) - nvars; % Number of data points minus one

    % Calculate reduced chi-squared
    if nu > 0
        fval = chi_squared / nu;
    else
        error('Not enough data points to calculate reduced chi-squared.');
    end
end

