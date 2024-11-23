function [Cp] = tima_specific_heat_model_DV1963(rho_dry,rho,Temp_K,material)
%% TIMA_VERHOEF_SPECIFIC_HEAT_MODEL_DV1963
%   function to calculates specific heat as a funtion of volumetric water content, temperature, and density
%
% Description
%   A Cp model that considers the change of Cp_solid with temperature 
%   and added the effects of solid and liquid in series (e.g., Verhoef 2004)
%
% Syntax
%   [Cp] = tima_verhoef_specific_heat_model(rho_dry,rho_H2O,rho,VWC,temperature)
%
% Inputs
%   rho_dry: [kg/m^3] dry soil density
%   rho: [kg/m^3] current soil density
%   VWC: [fraction by vol] current water content
%   Temp_K: [K] current temperature of soil
%
% Outputs
%   Cp: [%J/kgK] Specific heat capacity   
%
% Author
%    Ari Koeppel, 2021
%
% References
%   Waples 2004B: Somerton (1992) states that the
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
CpnT = @(T)8.95E-10.*T.^3-2.13E-6.*T.^2+0.00172.*T+0.716; %Waples2004a Temp in C
if strcmp(material,"basalt")
    Cp_solid = 898*CpnT(Temp_K-273.15)/CpnT(20); %basalt
elseif strcmp(material,"amorphous") %silica glass
    Cp_solid = 795*CpnT(Temp_K-273.15)/CpnT(20); %Volcanic tuff
elseif strcmp(material,"granite")
    Cp_solid = 1172*CpnT(Temp_K-273.15)/CpnT(20);
elseif strcmp(material,"sandstone")
    Cp_solid = 775*CpnT(Temp_K-273.15)/CpnT(20);
elseif strcmp(material,"clay")
    Cp_solid = 860*CpnT(Temp_K-273.15)/CpnT(20); %clay
elseif strcmp(material,"salt")
    Cp_solid = 880*CpnT(Temp_K-273.15)/CpnT(20); %Salt
elseif strcmp(material,"ice") %Waples 2004b
    Cp_solid = 7.8277.*(Temp_K-273.15)+2115;
else
    Cp_solid = -23.173+2.127.*Temp_K+1.5008e-2.*Temp_K.^2-7.3699e-5.*Temp_K.^3+9.6552E-8.*Temp_K.^4; %Cp of solids - polynomial fit for basalt mineral soils
end
Cp_H2O = 3.58435521E-06.*Temp_K.^4-0.00475626235.*Temp_K.^3+2.37137032.*Temp_K.^2-526.09.*Temp_K+47971; %Engineering Toolbox
Cp = (Cp_H2O.*(rho-rho_dry)+Cp_solid.*rho_dry)./rho; % combined by weight fraction
end