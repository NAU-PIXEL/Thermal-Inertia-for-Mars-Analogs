function [k_eff] = tima_conductivity_model_piqueux2009(k_dry_std,Soil_Temperature,T_std,VWC,theta_k,m,Soil_RH,material)
%% TIMA_CONDUCTIVITY_MODEL_piqueux2009
%   function to calculate the effective thermal conductivity of a particulate soil or sediment using the piqueux 2009 Method
%
% Description
%   [m = ~0.3 seems best] 0<m<1.33  Best for granular with no hydration regime, Realistic lower limit
%      The lu 2007 fitted values of m are 0.96 and 0.27 for the coarse-textured and fine textured soils, respectively.
%
% Syntax
%   [k_eff] = tima_conductivity_model_piqueux2009(k_dry_std,Soil_Temperature,T_std,VWC,theta_k,m,Soil_RH,material)
%
% Inputs
%   k_dry_std: [W/mK] stardard bulk dry thermal conductivty at T_std (vector)
%   Soil_Temperature: [K] Soil temperature (vector)
%   T_std: [K] Standard temperature; typically 300 (vector)
%   VWC: [fraction by volume] volumetric water content (vector)
%   k_H2O: [W/mK] Thermal conductivity of soil water (vector)
%   theta_k: [fraction by volume] conductivity soil moisture inflection point - in theory this should be similar to saturation value or porosity (vector)
%   m: [unitless] pore-network-connectivity-parameter (vector)
%   Soil_RH: [fraction] relative humidity of air within the soil (vector)
%   material: ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
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


% T_std = 300; % Temperature standard for thermal conductivity (K)

k_air = 0.024+7.78E-5.*(Soil_Temperature-273.15) - 3E-05.*exp(0.0591*(Soil_Temperature-273.15)).*Soil_RH; %Tsilingiris 2008 Fig 3 assuming soil is equilibrated with air
k_sc_mod = (1.840E-3+3.957E1*k_air+2.871E3*k_air^2+5.016E3*k_air^3)/(1+2.25E2*k_air+4.491E3*k_air^2+3.621E3*k_air^3);
k_cc_mod = (4.116E-5+1.108E1*k_air+3.229E2*k_air^2+5.366E2*k_air^3)/(1+7.407E1*k_air+6.248E2*k_air^2+2.095E2*k_air^3);
k_dry_mod = k_sc_mod+(k_cc_mod-k_sc_mod)*(theta_k-0.259)/(0.476-0.259);
if k_dry_mod<k_air
    k_dry_mod = k_air;
end


k_air_avg = 0.024+7.78E-5.*(T_std-273.15);
k_sc_mod_avg = (1.840E-3+3.957E1*k_air_avg+2.871E3*k_air_avg^2+5.016E3*k_air_avg^3)/(1+2.25E2*k_air_avg+4.491E3*k_air_avg^2+3.621E3*k_air_avg^3);
k_cc_mod_avg = (4.116E-5+1.108E1*k_air_avg+3.229E2*k_air_avg^2+5.366E2*k_air_avg^3)/(1+7.407E1*k_air_avg+6.248E2*k_air_avg^2+2.095E2*k_air_avg^3);
k_dry_mod_avg = k_sc_mod_avg+(k_cc_mod_avg-k_sc_mod_avg)*(theta_k-0.259)/(0.476-0.259);
if k_dry_mod_avg<k_air_avg
    k_dry_mod_avg = k_air_avg;
end

 k_dry = k_dry_avg*k_dry_mod/k_dry_mod_avg; % Since Piquex 2009a Fig 8 was modeled for a specific material type, k_dry_mod does not encessarily apply here, but the ratio of k_dry_mod/k_dry_mod_avg should be similar.
%     k_dry = k_dry_avg+3.66*(k_air-k_air_avg);

if strcmp(material,'basalt')
    k_solid = 1.18 + 474/(350+Soil_Temperature-273.15);%  % Piqueux and Christensen 2011/Clauser and Huenges [1995]
% Roughly = 2.2, Bristow, 2002
elseif strcmp(material,'amorphous')
    k_solid = 0.6924 + 0.0015*(Soil_Temperature-273.15); %Piqueux and Christensen 2011/Clauser and Huenges [1995]
elseif strcmp(material,'granite')
    k_solid = 2.0;
    % Roughly = 2.0 (granite), Bristow, 2002
elseif strcmp(material,'quartz')
    k_solid = 7.69; %(horai 1971)
elseif strcmp(material,'clay')
    k_solid = 2.9;
    % Roughly = 2.9 (clay), Bristow, 2002
elseif strcmp(material,'salt')
    k_solid = -2.11 + 2960/(350+Soil_Temperature-273.15); %Piqueux and Christensen 2011/Clauser and Huenges [1995]
elseif strcmp(material,'ice')
    k_solid = 2.18;
    % Roughly = 2.18, Bristow, 2002
else
    k_solid = 2.2;
end

k_H2O = -1.1e-5.*(Soil_Temperature - 273.15).^2 + 0.00234.*(Soil_Temperature - 273.15) + 0.552; %Bristow 2002 or 0.59 from Zhang and Wang 2017 or 0.6096 from Haigh 2012

Sr = VWC/theta_k;
k_wet = k_solid^(1-theta_k)*(k_H2O^theta_k);%geometric mean
if Sr>1, Sr = 1; end

k_eff = k_dry + exp(m*(1-(Sr)^(m-1.33)))*(k_wet-k_dry); %Lu 2007 [m = ~0.3 seems best] 0<m<1.33  Best for granular with no hydration regime, Realistic lower limit
end