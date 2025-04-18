function [q_latent,soil_RH] = tima_latent_heat_model_VDG1985(CE,theta_E,pressure_air_Pa,windspeed_horiz,RH,air_temp_K,T_ground_K,soil_VWC,CH)
%% TIMA_LATENT_HEAT_MODEL_VDG1985
%   function to calculate latent heat flux due to evaporation in a given
%   layer at a given timestep using a tunable aerodynamic method from Van de Griend et al., 1985
%   **This script is incomplete (Jan, 2025)**
%
% Syntax
%   [q_latent,Soil_RH] = tima_latent_heat_model_VDG1985(1000,0.2,1000,5,50,300,290,0.1,250)
%
% Inputs
%   air_temp_K: [K] near surface air temperature, typically at 3 m AGL.
%           (scalar)
%   CE: [Unitless] resistance to latent heat flux coefficient, similar to the
%       aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2)
%       (scalar)
%   CH: [Unitless] resistance to sensible heat flux coefficient,
%           similar to the aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2)
%           1/CH should be 0.0028-0.0075 (or CH~100-400) for smooth to
%           roughly open soils on Davenport Scale, CH is larger (more
%           resistence) with rougher topography or larger vegetation (scalar)
%   pressure_air_Pa: [Pa] station pressure, typically at 3m AGL. (scalar)
%   RH: [0-1] array of near surface relative humidity values,
%           typically at 3m AGL. (scalar)
%   soil_VWC: [0-1] soil volumetric water content at depth of interest (scalar)
%   T_ground_K: [K] soil temperature at depth of interest (scalar)
%   theta_E: [0-1, fraction by volume] latent heat soil moisture inflection
%           point (scalar)
%   windspeed_horiz: [m/s] Near surface horizontal wind speed,
%           typically at 3m AGL. (scalar)
%
% Outputs
%   q_latent: [W/m^2] latent heat flux from evaporation (scalar)
%   soil_RH: [0-1] Theoretical "relative humidity" within soil (scalar)
%
% Author
%    Ari Koeppel, 2021
%
% Sources
%   Campbell 1985
%   Price 1985
%   van de Griend et al 1985
%   http://cires1.colorado.edu/~voemel/vp.html BEST REVIEW OF SATURATION VAPOR PRESSURES
%   Evett in Warrick 2002
%   Sonntag 1994
%   Mahfouf & Noilhan, 1991
%   Tran 2016
%   van de Griend and Owe 1994
%   Camillo and Gurney 1986
%   Fen Shu 1982
%   Daamen and Simmonds 1996

    if windspeed_horiz < 1
        windspeed_horiz = 1; %m/s, lower limit on measureable wind speed 
    end
    L_H2O = (2.501-2.37e-3*(T_ground_K-273.15))*10^6;%J/kg, latent heat of evaporation  - Evett in Warrick 2002
    R = 8.31446; %J/(K·mol), gas constant
    Md = 0.0289652; % kg/mol, Dry air
    Mv = 0.018016; % kg/mol, H2O
    Cp_air = 1013;%J⋅kg−1⋅K−1 Evett in Warrick 2002
    Psat_air = 100*exp(-6096.9385/air_temp_K+16.635794-2.711193e-2*air_temp_K+1.673952e-5*air_temp_K^2+2.433502*log(air_temp_K)); %Pa, Sonntag 1994 see http://cires1.colorado.edu/~voemel/vp.html
    Pv = Psat_air*RH; %Pa, Part Press vapor in air
    if isempty(pressure_air_Pa)
        rho_air = 1.291-0.00418*(air_temp_K-273.15); %kg/m^3, Evett in Warrick 2002
    else
        Pd = pressure_air_Pa - Pv; %Pa, Part Press dry air
        rho_air = (Pd*Md+Pv*Mv)/(R*air_temp_K); %kg/m^3, accounting for humid air
    end
    Psat_surf = 100*exp(-6096.9385/T_ground_K+16.635794-2.711193e-2*T_ground_K+1.673952e-5*T_ground_K^2+2.433502*log(T_ground_K)); %Pa, Sonntag 1994 see http://cires1.colorado.edu/~voemel/vp.html
    q_air = Mv/Md*Pv/(pressure_air_Pa-(1-Mv/Md)*Pv); %kgh2o/kg dry_air, Specific Humidity Air
    q_sat_surf = Mv/Md*Psat_surf/(pressure_air_Pa-(1-Mv/Md)*Psat_surf); %kgh2o/kg dry_air, Specific Humidity Surf (saturation value at surf temp)
    if soil_VWC <= theta_E % theta_E is max amount a soil will hold if exposed to a constant mist, Mahfouf & Noilhan, 1991, Tran 2016 calls VWC_Sat the evaporation-rate reduction point
        Rs = (theta_E-soil_VWC)*CE;
    else
        Rs = 0;
    end
    if q_sat_surf<q_air
       soil_RH = RH; % Theoretical "relative humidity" within soil, Pct by vol, Lee & Pielke 1992 + Tran, 2016
    else
       soil_RH = 1; % saturated
    end
    Ra = (rho_air*Cp_air)/(CH*windspeed_horiz);
    %
%     Rs = 3e10*(VWC_Sat-Soil_VWC)^16.6; %Daamen and Simmonds 1996
%     Rs = -885+4140*(VWC_Sat-Soil_VWC);  %Soil_VWC in top 0.5 cm Camillo and Gurney 1986
%     Rs = 3.5*(VWC_Sat/Soil_VWC)^2.3+33.5;%Fen Shu 1982 top 5 mm
%     Rs = 10*exp(0.3565*(VWC_Sat-Soil_VWC));  %Soil_VWC in top 1 cm van de Griend and Owe 1994 + Tran 2016
    Evap = rho_air*(q_sat_surf-q_air)/(Ra + Rs);
    q_latent = -L_H2O*Evap; %W/m^2, Latent heat flux

end
