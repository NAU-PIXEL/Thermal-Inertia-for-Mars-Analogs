function [logfval] = tima_logPfun_dchi2v_dt(observed,modeled,error)
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
    dobs = gradient(observed);
    dmod = gradient(modeled);

    % Remove the last element to match sizes (ddata1 and ddata2 will be 1 shorter)
    dobs = movmean(dobs(1:end-1),5);
    dmod = movmean(dmod(1:end-1),5);
    unc = movmean(error(1:end-1),5);

    logfval = -0.5*sum((dobs-dmod).^2./(unc).^2);% a cell of function handles returning the log probality of a each outcom
end

