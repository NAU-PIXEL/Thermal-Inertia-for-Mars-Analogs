function [logp] = tima_ln_prior(FitParams,k_upper_min,k_upper_max,k_lower_min,k_lower_max,depth_min,depth_max,m_min,m_max,CH_min,CH_max,CE_min,CE_max,theta_k_min,theta_k_max,theta_E_min,theta_E_max)
%***************
% TIMA_LN_PRIOR
%    Set up resonable ranges for each test variable. 
% 
% INPUTS
%   FitParams: parameters for fitting
%   min and max for
%   k (upper and lower) = thermal conductivity (W/mK); k_air = 0.024 from Tsilingiris 2008; k_solid_basalt = 2.2 from Bristow, 2002; k_solid_granite = 3.7 from Piqueux 2011;
%   depth = subsurface depth transition between k_upper and k_lower (m)
%   m = pore-network-connectivity (unitless)
%   CH = sensible heat resistence coefficent (unitless)
%   CE = latent heat resistence coefficient (unitless)
%   theta_k = soil moisture inflection point for thermal conductivity (% by volume)
%   theta_E = soil moisture inflection point for latent heat (% by volume)
% 
% OUTPUT 
%   logp = If the simulation is within the ranges, the log-probability is 0. Otherwise, it's -infinity.
%
% AUTHOR
%   Ari Koeppel, 2023

    p = inputParser;
    p.addRequired('theta',@(x) all(isvector(x) & isnumeric(x)));
    p.addOptional('k_upper_min',0.024, @(x) all(isvector(x) & isnumeric(x) & x>0)); 
    p.addOptional('k_upper_max',3.7, @(x) all(isvector(x) & isnumeric(x) & x<10 & x>0)); 
    k_upper_min = p.Results.k_upper_min(:);
    k_upper_max = p.Results.k_upper_max(:);
    p.addOptional('k_lower_min',0.024, @(x) all(isvector(x) & isnumeric(x) & x>0)); 
    p.addOptional('k_lower_max',3.7, @(x) all(isvector(x) & isnumeric(x) & x<10 & x>0)); 
    k_lower_min = p.Results.k_lower_min(:);
    k_lower_max = p.Results.k_lower_max(:);
    p.addOptional('depth_min',0, @(x) all(isvector(x) & isnumeric(x) & x>=0));
    p.addOptional('depth_max',1.3, @(x) all(isvector(x) & isnumeric(x) & x<=1.3 & x>0));
    depth_min = p.Results.depth_min(:);
    depth_max = p.Results.depth_max(:);
    p.addOptional('m_min',0, @(x) all(isvector(x) & isnumeric(x) & x>=0));
    p.addOptional('m_max',1.3, @(x) all(isvector(x) & isnumeric(x) & x<=1.3 & x>0));
    m_min = p.Results.m_min(:);
    m_max = p.Results.m_max(:);
    p.addOptional('CH_min',1, @(x) all(isvector(x) & isnumeric(x) & x>=0));
    p.addOptional('CH_max',10000, @(x) all(isvector(x) & isnumeric(x) & x>0));
    CH_min = p.Results.CH_min(:);
    CH_max = p.Results.CH_max(:);
    p.addOptional('CE_min',1, @(x) all(isvector(x) & isnumeric(x) & x>=0));
    p.addOptional('CE_max',10000, @(x) all(isvector(x) & isnumeric(x) & x>0));
    CE_min = p.Results.CE_min(:);
    CE_max = p.Results.CE_max(:);
    p.addOptional('theta_k_min',0.05, @(x) all(isvector(x) & isnumeric(x) & x>=0 & x<1));
    p.addOptional('theta_k_max',0.9, @(x) all(isvector(x) & isnumeric(x) & x<=1 & x>0));
    theta_k_min = p.Results.theta_k_min(:);
    theta_k_max = p.Results.theta_k_max(:);
    p.addOptional('theta_E_min',0.01, @(x) all(isvector(x) & isnumeric(x) & x>=0 & x<1));
    p.addOptional('theta_E_max',0.75, @(x) all(isvector(x) & isnumeric(x) & x<=1 & x>0));
    theta_E_min = p.Results.theta_E_min(:);
    theta_E_max = p.Results.theta_E_max(:);

    if (FitParams(1)>k_upper_min) && (FitParams(1)<k_upper_max)...
            && (FitParams(2)>k_lower_min) && (FitParams(2)<k_lower_max)...
            && (FitParams(3)>depth_min) && (FitParams(3)<depth_max)...
            && (FitParams(4)>m_min) && (FitParams(4)<m_max)...
            && (FitParams(5)>CH_min) && (FitParams(5)<CH_max) && (FitParams(6)>CE_min)...
            && (FitParams(6)<CE_max) && (FitParams(7)>=theta_k_min) && (FitParams(7)<theta_k_max)...
            && (FitParams(8)>theta_E_min) && (FitParams(8)<theta_E_max) 
        logp = 0;
    else
        logp = -Inf;
    end
end