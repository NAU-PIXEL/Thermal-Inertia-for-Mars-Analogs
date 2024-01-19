function [Subsurface_Temperatures_C] = tima_initialize(kay,rho_dry,m,VWC_Sat,T_Deep,Soil_Temp_C,dt,B,Dug_VWC,VWC_depth_indices,RH)
%***************
% TIMA_INITIALIZE
%   This function uses observational data and assigned thermophysical properties 
%   to estimate the emitting surface temperature through time.

% Syntax
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location)
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location, Pressure)
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location, Pressure, Temperature)
%   [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location, 'temperature', Temperature)

% Description
%   Simple version of Heat tansfer model used to estimate realistic subsurface 
%   temperatures at start of model using the measured surface temp for 1 day on repeat as the forcing.
%
% Input Parameters
%   Observational data
%       dt,
%       Air_Temp_C
%       GHI
%       R_Short_Lower
%       R_Long_Upper
%       WindSpeed_horiz
%       Dug_VWC
%       VWC_depth_indices
%       RH
%       rho_dry
%       emissivity
%       albedo
%       pressure_air
%    Boundary conditions
%       ST1
%       Initial_Temps
%       T_Deep
%       Cp_Deep
%       B
%       dt
%    Fitting Parameters (each may be fit or assigned)
%       kay_upper
%       kay_lower
%       Depth_transition
%       m
%       Wind_Coeff
%       Evap_Coeff
%       theta_k
%       theta_E
%
% Outputs:
%   T_Surf_C = Surface temperature time series (deg C)
%
% Author
%    Ari Koeppel -- Copyright 2023
%
% Sources
%   Optimization tool: https://www.mathworks.com/help/gads/table-for-choosing-a-solver.html
%   Subsurface heat flux: Kieffer et al., 2013
%   Irradiance Calculation: https://github.com/sandialabs/MATLAB_PV_LIB
%   Sensible heat flux:
%   Latent heat flux: Daamen and Simmonds 1996 + Mahfouf and Noilhan (1991)+ Kondo1990
%   Thermal conductivity mixing:
%   
% See also 
%   TIMA_HEAT_TRANSFER TIMA_INITIALIZE TIMA_LATENT_HEAT_MODEL TIMA_LN_PRIOR TIMA_SENSIBLE_HEAT_MODEL TIMA_GWMCMC TIMA_COMBINE_ROWS
    % ***************    
    
    % ***************  
    % # of days to run equilib model for: More = closer to equilib, fewer = faster
    NDAYS = 20; %run the inputs from the first day 20 times
    % ***************  
    
    % ***************
    % Inputs and constants:
    Soil_Temp_K = Soil_Temp_C+273.15;
    Day_Dur = 1440/(dt/60); %#of mins/day div by min per interval
    Depth = sum(B);
    NLAY = length(B);
    k = kay.*ones(1,NLAY);
    k(NLAY) = kay;
    rho_H2O = 997;
    k_H2O = 0.6096; %from Haigh 2012
    % ***************

    % ***************
    % Initialize temperatures of subsurface layers linearly at start
    Subsurface_Temperatures_C = NaN(NDAYS,Day_Dur,NLAY); %Set up Matrix of Temperatures for each day, each minute, and each layer
    Subsurface_Temperatures_C(:,:,NLAY) = T_Deep; %Constant lower boundary temp
    for ind = 2:NLAY-1
        Subsurface_Temperatures_C(1,1,ind) = Soil_Temp_K(1)+(T_Deep-Soil_Temp_K(1))*(sum(B(1:ind))/Depth);
        k(ind) = Conductivity_model(kay,Subsurface_Temperatures_C(1,1,ind),T_Deep,Dug_VWC(1,VWC_depth_indices(end)),k_H2O,VWC_Sat,m,rho_dry,RH(1)); %Temp and H2O dep
    end
    % ***************
    
    % ***************
    % Cycle through day 1 NDAYS times using measured surf temp and VWCs to estimate subsurface forcing
    for day = 1:NDAYS %Day loop
        Subsurface_Temperatures_C(day,:,1) = Soil_Temp_K(1:Day_Dur); %Layer 1 temp = observed
        for t = 1:Day_Dur %Minute loop
            k(1) = 1/(Dug_VWC(t,VWC_depth_indices(1))/k_H2O+(kay*sqrt(Subsurface_Temperatures_C(day,t,1)/T_Deep))^-1);
            for i = 2:NLAY-1 %layer loop
                if t == 1 && day > 1 %First step in day (requires pulling from previous day)
                    rho = rho_dry + rho_H2O*Dug_VWC(Day_Dur,VWC_depth_indices(i)); %H2O dep
                    Cp = Specific_Heat_Model(rho_dry,rho_H2O,rho,Dug_VWC(Day_Dur,VWC_depth_indices(i)),Subsurface_Temperatures_C(day-1,end,i));
                    k(i) = Conductivity_model(kay,Subsurface_Temperatures_C(day-1,end,i),T_Deep,Dug_VWC(Day_Dur,VWC_depth_indices(i)),k_H2O,VWC_Sat,m,rho_dry,RH(Day_Dur));
                    %Subsurface Multip Factors (Kieffer, 2013)
                    F1 = 2*dt*k(i)/((Cp*rho*B(i)^2)*(1+B(i+1)/B(i)*k(i)/k(i+1)));
                    F3 = (1+(B(i+1)*k(i))/(B(i)*k(i+1)))/(1+(B(i-1)/B(i)*k(i)/k(i-1)));
                    F2 = -(1 + F3);
                    %Temperature response
                    dT = F1*(Subsurface_Temperatures_C(day-1,end,i+1)+F2*Subsurface_Temperatures_C(day-1,end,i)+F3*Subsurface_Temperatures_C(day-1,end,i-1));
                    Subsurface_Temperatures_C(day,t,i)  = Subsurface_Temperatures_C(day-1,end,i) + dT;
                elseif t >= 2 %All other steps
                    rho = rho_dry + rho_H2O*Dug_VWC(t-1,VWC_depth_indices(i)); %H2O dep
                    Cp = Specific_Heat_Model(rho_dry,rho_H2O,rho,Dug_VWC(t-1,VWC_depth_indices(i)),Subsurface_Temperatures_C(day,t-1,i));
                    k(i) = Conductivity_model(kay,Subsurface_Temperatures_C(day,t-1,i),T_Deep,Dug_VWC(t-1,VWC_depth_indices(i)),k_H2O,VWC_Sat,m,rho_dry,RH(t-1));
                    %Subsurface Multip Factors (Kieffer, 2013)
                    F1 = 2*dt*k(i)/((Cp*rho*B(i)^2)*(1+B(i+1)/B(i)*k(i)/k(i+1)));
                    F3 = (1+(B(i+1)*k(i))/(B(i)*k(i+1)))/(1+(B(i-1)/B(i)*k(i)/k(i-1)));
                    F2 = -(1 + F3);
                    %Temperature response
                    dT = F1*(Subsurface_Temperatures_C(day,t-1,i+1)+F2*Subsurface_Temperatures_C(day,t-1,i)+F3*Subsurface_Temperatures_C(day,t-1,i-1));
                    Subsurface_Temperatures_C(day,t,i) = Subsurface_Temperatures_C(day,t-1,i)+dT;
                end                 
            end
        end
    end
end