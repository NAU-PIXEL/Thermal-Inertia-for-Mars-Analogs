function [k_eff] = tima_conductivity_model_chen2008(VWC,theta_k,m)
%% TIMA_CONDUCTIVITY_MODEL_chen2008
%   function to calculate the effective thermal conductivity of a particulate soil or sediment using the Chen 2008 Method
%
% Description
%   0<m<1
%
% Syntax
%   [k_eff] = tima_conductivity_model_chen2008(k_dry_std,Soil_Temperature,T_std,VWC,k_H2O,theta_k,m,Soil_RH,material)
%
% Inputs
%   VWC: [fraction by volume] volumetric water content (vector)
%   theta_k: [fraction by volume] conductivity soil moisture inflection point - in theory this should be similar to saturation value or porosity (vector)
%   m: [unitless] pore-network-connectivity-parameter (vector)
%
% Outputs
%   k_eff: [W/mK] effective thermal conductivity (vector)
%
% Author
%    Ari Koeppel, 2021
%
% Sources
%   Zhang and Wang 2017, Dong 2015 reviews & Lu and Dong 2015
%   Tsilingiris 2008 + Piqueux 2009a
%   Piqueux and Christensen 2011
%   Clauser and Huenges [1995]
%   Bristow, 2002: Basalt: 2.2, Granite: 2.0, Quartz: 8.8, Clay: 2.9, organics: 0.25, Ice @ 0C: 2.18
%   Water = 0.552+2.34E-3*Soil_Temperature-1.1E-5*Soil_Temperature^2
%   Air: 0.0237 + 0.000064*Soil_Temperature
%   Campbell 1994: Dry_air = 0.024+7.73E-5.*(Soil_Temperature-273.15) - 2.6E-8*(Soil_Temperature-273.15)^2
%   Horai 1971
%   Pielke 2002; k_eff Concrete: 4.6, Rock: 2.93, Ice: 2.51, Snow: 0.08-1.67, Stable air: 0.02-0.03, water; 0.57-0.63
%   Clay-dry: 0.25, clay-10%h2o: 0.63, clay-20%h2o: 1.12, clay-30%h2o: 1.33, clay-40%h2o: 1.58, --porosity 40%
%   sand-dry: 0.30, sand-10%h2o: 1.05, sand-20%h2o: 1.95, sand-30%h2o: 2.16, sand-40%h2o: 2.20, --porosity 40%
%   peat-dry: 0.06, peat-10%h2o: 0.10, peat-20%h2o: 0.29, peat-30%h2o: 0.43, peat-40%h2o: 0.50, --porosity 80%
%   rooty-soil: 0.11;
%   Lee & Pielke 1992: Field capacity- Sand: 0.135, Loam: 0.255, Clay: 0.367, Peat: 0.535
%   Saturation- Sand: 0.395, Loam: 0.451, Clay: 0.482, Peat: 0.863


Sr = VWC/theta_k; %degree of field capacity

if Sr>1, Sr = 1; end

%**********K Model************
k_eff = 7.5^(1-porosity)*0.61^(theta_k)*((1-m)*Sr+m)^(0.78*theta_k); %Chen 2008
end