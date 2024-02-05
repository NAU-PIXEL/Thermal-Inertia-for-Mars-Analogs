function [Cp] = tima_specific_heat_model_hillel(rho_dry,rho,VWC)
% TIMA_HILLEL_SPECIFIC_HEAT_MODEL
%   function to calculates specific heat as a funtion of volumetric water content
%
% Description
%   A Cp model that considers water content and yields similar values as Verhoef model (considering temperature) while saving computation time.
%
% Syntax
%   [Cp] = tima_specific_heat_model_hillel(rho_dry,rho,VWC)
%
% Inputs
%   rho_dry: [kg/m^3] dry soil density, for ice set to 950
%   rho: [kg/m^3] current soil density, for ice set to 950
%   VWC: [fraction by vol] current water content, for ice set to 0.3
%
% Outputs
%   Cp: [%J/kgK] Specific heat capacity   
%
% Author
%    Ari Koeppel, 2021
%
% Sources
%   Hillel 1980 + Evett in Warrick 2002
%   Pielke 2002; Concrete: 879, Rock: 753, Ice: 2100, Snow: 2093, Stable air: 1005, water; 4186
%   Clay-dry: 890, clay-10%h2o: 1005, clay-20%h2o: 1172, clay-30%h2o: 1340, clay-40%h2o: 1550,  --porosity 40%
%   sand-dry: 800, sand-10%h2o: 1088, sand-20%h2o: 1256, sand-30%h2o: 1423, sand-40%h2o: 1480,  --porosity 40%
%   peat-dry: 1920, peat-10%h2o: 2302, peat-20%h2o: 3098, peat-30%h2o: 3433, peat-40%h2o: 3650, --porosity 80%
%   rooty-soil: 1256;
%   Abu-Hamdeh -- clay: 1170-2250, sand: 830-1670

Cp = (2E6.*rho_dry./1000./2.65+4.2E6.*VWC)./rho;
end