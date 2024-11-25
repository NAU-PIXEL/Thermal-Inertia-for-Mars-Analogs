function [Psat] = tima_psaturationpa_Sonntag1994(TempK)
%% TIMA_PSATURATIONPA_SONNTAG1994
%   function to saturation vapor pressure at a given temperature
%
% Syntax
%   [Psat] = tima_PSaturationPa_Sonntag1994(TempK)
%
% Inputs
%   TempK: temperature (K)
%
% Outputs
%   Psat: Saturation vapor pressure (Pa)
%
% Sources:
%   Sonntag 1994
%   http://cires1.colorado.edu/~voemel/vp.html

Psat = 100.*exp(-6096.9385./TempK+16.635794-2.711193E-2.*TempK+1.673952E-5*TempK.^2+2.433502.*log(TempK)); %Pa
end