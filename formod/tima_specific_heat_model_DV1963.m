function [Cp] = tima_specific_heat_model_DV1963(rho_dry,rho,Temp_K,material)
%% TIMA_SPECIFIC_HEAT_MODEL_DV1963
%   function to calculate specific heat as a funtion of water content,
%   temperature, mineralogy and density, in J/kgK
%
% Description
%   A Cp model that considers the change of Cp_solid with temperature 
%   and added the effects of solid and liquid in series (e.g., Verhoef 2004)
%
% Syntax
%   [Cp] = tima_verhoef_specific_heat_model(1300,1600,300,"basalt")
%
% Inputs
%   rho_dry: [kg/m^3] dry soil density (scalar)
%   rho: [kg/m^3] current soil density (scalar)
%   Temp_K: [K] current temperature of soil (scalar)
%
% Outputs
%   Cp: [J/kgK] Specific heat capacity (scalar)
%
% Author
%    Ari Koeppel, 2021
%
% References
%   Waples and Waples 2004B: Somerton (1992) states that the
%       contribution from gas can be ignored without loss of
%       accuracy, but it is easy to include that of gas if desired.
%   CS HFP01 Soil Heat Flux Plate Manual + Verhoef 2004 -- Cp dry mineral
%   soil = 840.
%   Wang 2016
%   Abu-Hamdeh 2003
%   de Vries 1963
%   Pielke 2002; Concrete: 879, Rock: 753, Ice: 2100, Snow: 2093, Stable air: 1005, water; 4186
%   Clay-dry: 890, clay-10%h2o: 1005, clay-20%h2o: 1172, clay-30%h2o: 1340, clay-40%h2o: 1550,  --porosity 40%
%   sand-dry: 800, sand-10%h2o: 1088, sand-20%h2o: 1256, sand-30%h2o: 1423, sand-40%h2o: 1480,  --porosity 40%
%   peat-dry: 1920, peat-10%h2o: 2302, peat-20%h2o: 3098, peat-30%h2o: 3433, peat-40%h2o: 3650, --porosity 80%
%   rooty-soil: 1256;
%   Bristow 2002
%   Abu-Hamdeh -- clay: 1170-2250, sand: 830-1670
%   Solid Rock Eq from Robertson, E.C. and Hemingway, B.S., 1995. Estimating heat capacity and heat content of rocks (No. 95-622). US Geological Survey.

CpnT = @(T)8.95E-10.*T.^3-2.13E-6.*T.^2+0.00172.*T+0.716; %normalized heat capacity (Cpn) of a
    %mineral or a nonporous rock at any temperature, Waples2004a Temp in C
%% Equations below use the ratio of the normalized heat capacity at a given temperature to a known value at 20 C (except ice), J/kgK
if strcmp(material,"basalt")
    Cp_solid = 898*CpnT(Temp_K-273.15)/CpnT(20);
elseif strcmp(material,"amorphous") %silica glass/volcanic tuff
    Cp_solid = 795*CpnT(Temp_K-273.15)/CpnT(20);
elseif strcmp(material,"granite")
    Cp_solid = 1172*CpnT(Temp_K-273.15)/CpnT(20);
elseif strcmp(material,"sandstone")
    Cp_solid = 775*CpnT(Temp_K-273.15)/CpnT(20);
elseif strcmp(material,"clay")
    Cp_solid = 860*CpnT(Temp_K-273.15)/CpnT(20);
elseif strcmp(material,"salt")
    Cp_solid = 880*CpnT(Temp_K-273.15)/CpnT(20); 
elseif strcmp(material,"ice") %Waples 2004b
    Cp_solid = 7.8277.*(Temp_K-273.15)+2115;
else
    Cp_solid = -23.173+2.127.*Temp_K+1.5008e-2.*Temp_K.^2-7.3699e-5.*Temp_K.^3+9.6552E-8.*Temp_K.^4; %Polynomial fit for typical mineral soils
end
Cp_H2O = 3.58435521E-06.*Temp_K.^4-0.00475626235.*Temp_K.^3+2.37137032.*Temp_K.^2-526.09.*Temp_K+47971; %Specific Heat of liquid water. J/kgK
Cp = (Cp_H2O.*(rho-rho_dry)+Cp_solid.*rho_dry)./rho; %Equation from De Vries 1963. J/kgK
end