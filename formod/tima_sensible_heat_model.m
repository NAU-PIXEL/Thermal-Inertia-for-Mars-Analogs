function [q_sensible,SHCoeff] = tima_sensible_heat_model(CH,windspeed_horiz,air_temp_K,T_surf_K)
%% TIMA_SENSIBLE_HEAT_MODEL
%   function to calculate sensible heat flux at the surface due to forced convection at a given timestep using a tunable aerodynamic method
%
% Syntax
%   [q_sensible] = tima_sensible_heat_model(CH,windspeed_horiz,air_temp_K,T_surf_K)
%
% Inputs
%   CH: resistance to sensible heat flux coefficient, similar to the aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2) (Unitless)
%       1/CH should be 0.0028-0.0075 (or CH~100-400) for smooth to roughly open soils on Davenport Scale, CH is larger (more resistence) with rougher topography or larger vegetation
%   windspeed_horiz: windspeed from tower height associated with derived CE (m/s)
%   air_temp_K: air temperature from tower height associated with derived CE (K)
%   T_surd_K: ground surface temperature at depth of interest (K)
%
% Outputs
%   q_sensible: sensible heat flux from evaporation (W/m^2)
%
% Author
%    Ari Koeppel, 2021
%
% Sources
%   Cellier 1996
%   12th Amer. Meteorol. Soc. Symposium on Applied Climatology, 2000, pp. 96–99.

    if windspeed_horiz < 1
        windspeed_horiz = 1; %m/s, lower limit on measureable wind speed
    end
    rho_air = 1.291-0.00418*(air_temp_K-273.15); %kg/m^3, Evett in Warrick 2002
    Cp_air = 1013;%J/(kg⋅K), Evett in Warrick 2002
    SHCoeff=-rho_air*Cp_air*windspeed_horiz/CH;
    q_sensible = SHCoeff*(T_surf_K-air_temp_K); %W/m^2
end