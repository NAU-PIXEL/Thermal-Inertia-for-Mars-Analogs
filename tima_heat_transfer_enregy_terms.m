function [T_Surf_C] = tima_heat_transfer(k_dry_300_upper,k_dry_300_lower,depth_transition,...
    m,CH,CE,theta_k,theta_E,rho_dry_upper,rho_dry_lower,dt,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,VWC_depth_indices,RH,emissivity,...
    pressure_air_pa,albedo,aspect_cwfromS,solar_azimuth_cwfromS,solar_zenith_apparent,f_diff,MappingMode,material)
%***************
% TIMA_HEAT_TRANSFER
%   This function uses observational data and assigned thermophysical properties 
%   to estimate the emitting surface temperature through time.

% Syntax
%   [T_Surf_C] = tima_heat_transfer(k_dry_300_upper,k_dry_300_lower,depth_transition,...
    % m,CH,CE,theta_k,theta_E,rho_dry_upper,rho_dry_lower,dt,air_temp_C,GHI,r_Long_Upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,VWC_depth_indices,RH,emissivity,...
    % pressure_air_pa)
%   [T_Surf_C] = tima_heat_transfer(k_dry_300_upper,k_dry_300_lower,depth_transition,...
    % m,CH,CE,theta_k,theta_E,rho_dry_upper,rho_dry_lower,dt,air_temp_C,GHI,r_Long_Upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,VWC_depth_indices,RH,emissivity,...
    % pressure_air_pa,r_short_lower,albedo,aspect_cwfromS,solar_azimuth_cwfromS,solar_zenith_apparent,DF,MappingMode)

% Description
%   This model uses a surrogate optimization solver with 1-7 paramters to fit IR-image-observed
%   surface temperatures using meteorological data and surface parameters derived from
%   aerial/satellite imagery. Fitting Parameters can include any or all of: top layer thermal
%   conductivity at 300K, bottome layer thermal conductivity at 300K, Depth of thermal conductivity
%   transition, pore-network-connectivity, Surf. ex. coef. (sensible), Surf. ex. coef. (latent), Soil Moist. Inflection (conductivity) (%), Soil Moist. Infl. (latent heat) (%)
%
% Input Parameters
%   Mapping Mode:
%   Observational data
%       dt
%       Air_Temp_C
%       R_Short_Upper: must be positive
%       R_Short_Lower: Must be positive
%       R_Long_Upper
%       WindSpeed_horiz
%       Dug_VWC
%       VWC_depth_indices
%       RH
%       rho_dry
%       emissivity
%       albedo
%       pressure_air
%       f_diff
%       material
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
p = inputParser;
% p.addRequired('Time',@isstruct);
% p.addRequired('Location',@isstruct);
p.addOptional('MappingMode',0);
p.addOptional('slope_angle',0);
p.addOptional('albedo',0);
p.addOptional('aspect',0);
p.addOptional('solar_azimuth_cwfromS',[]);
p.addOptional('solar_zenith_apparent',[]);
p.addOptional('DF',[]);
p.addOptional('k_lower',k_dry_300_upper);
p.addOptional('rho_dry_lower',rho_dry_lower);
p.addOptional('depth_transition',sum(layer_size));
p.parse(varargin{:});

r_short_lower = p.Results.r_short_lower;
MappingMode = p.Results.MappingMode;
slope_angle = p.Results.slope_angle;
albedo = p.Results.albedo;
aspect_cwfromS = p.Results.aspect_cwfromS;
rho_dry_lower = p.Results.rho_dry_lower;
k_dry_300_lower = p.Results.k_dry_300_lower;
depth_transition = p.Results.depth_transition;



% ***************

% ***************
% Inputs and constants:
air_temp_K = air_temp_C+273.15;
% k_H2O = 0.6096; %from Haigh 2012
rho_H2O = 997; %kg/m^3, does not vary much in T range so assume constant density
sigma = 5.670374419e-8; % Stefan-Boltzmann Const
omega = (1+cosd(slope_angle))/2; % visible sky fraction ISOtropic sky model
T_std = 300; %Standard temperature K
% ***************
% Initial and boundary conditions setup + array initiation:
NLAY = length(layer_size);
T = NaN(length(air_temp_C),NLAY); %Matrix of Temperatures for each layer and time point, T(1,t) is surface temp array that will be fit to obs data
T(1,:) = initial_temps;
k = k_dry_300_upper.*ones(1,NLAY);
kay_top = k_dry_300_upper;
kay_bottom = k_dry_300_lower;
kay_upper = k_dry_300_upper;
kay_lower = k_dry_300_lower;
[~,D_i] = min(abs(depth_transition-cumsum(layer_size))); %index of depth for k transition
T(:,D_i:NLAY) = T_deep; %Constant lower boundary temp
    
if MappingMode == 0 %Tower mode
    for t = 2:length(air_temp_C)
       
        %*********LATENT HEAT***********
        [q_evap_1,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,1),dug_VWC(t-1,1));
        %*******************************

        %*********THERMAL CONDUCTIVITY***********
        k(1) = tima_conductivity_model_lu2007(kay_top,T(t-1,1),T_std,dug_VWC(t-1,1),theta_k,m,Soil_RH, material);%Top
        k(NLAY) = kay_bottom;%Bottom
        %****************************************
        
        %*********Specific heat***********
        Cp = tima_specific_heat_model_hillel(rho_dry,rho,dug_VWC(t-1,1));
        %****************************************

        %*********SENSIBLE HEAT***********
        q_conv = tima_sensible_heat_model(CH,windspeed_horiz(t-1),air_temp_K(t-1),T(t-1,1));
        %*******************************
              
               
        %*********Radiative Flux***********
        phi = solar_azimuth_cwfromS(t-1);
        Z = solar_zenith_apparent(t-1);
        Incidence = cosd(Z)*cosd(slope_angle)+sind(Z)*sind(slope_angle)*cosd(phi-aspect_cwfromS);Incidence(Incidence>1) = 1; Incidence(Incidence<-1) = -1;
        q_rad_full = (1-albedo)*Incidence*(1-f_diff)*r_short_upper(t-1)/cosd(Z)+(1-albedo)*omega*f_diff*r_short_upper(t-1)+(1-albedo)*(1-omega)*r_short_lower(t-1)+omega*emissivity*r_long_upper(t-1)-omega*emissivity*sigma*T(t-1,1)^4; %net heat flux entering surface
        Lit_Fraction = shadow_data(shadow_time_ind(t-1)); %Assumes shadow = 0, and unshadowed = 1;
        q_rad = Lit_Fraction*q_rad_full+(1-Lit_Fraction)*((1-albedo)*omega*f_diff*r_short_upper(t-1)+(1-albedo)*(1-omega)*r_short_lower(t-1)+omega*emissivity*r_long_upper(t-1)-omega*emissivity*sigma*T(t-1,1)^4); %net heat flux entering surface
        %*******************************
        
        %*********Ground Flux***********
        %ref: Kieffer, 2013
        q_G = -k(1)*(T(t-1,1)-T(t-1,2))/layer_size(1); %heat flux from conducting with lower layer
        %*******************************

        %*********Combine Heat transfer elements***********
        W = q_G+q_conv+q_rad+q_evap_1; %Heat Flux, W/m^2
        rho = rho_dry + rho_H2O*dug_VWC(t-1,1); %Dry density + water content
        dT = dt/(Cp*rho*layer_size(1))*(W); %Temperature change for given thickness of material with known volumetric Cp
        T(t,1) = T(t-1,1) + dT; %new T at surface
        %*******************************

        %*********Subsurface Flux (from KRC, Kieffer, 2013)***********
        for i = 2:NLAY-1 %layer loop
%           k_H2O = 1.1e-05*(T(t-1,1)- 273.15).^2 + 0.00234*(T(t-1,1)- 273.15) + 0.552; %Bristow 2002, temp dep
            if i < D_i
                if sum(layer_size(1:i)) <= 0.05 %Wang 2016 inflection point
                    [q_evap_i,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,i),dug_VWC(t-1,VWC_depth_indices(i)));
                else
                    [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,i),dug_VWC(t-1,VWC_depth_indices(i)));
                    q_evap_i = 0; %Evap_Coeff
                end
                k(i) = tima_conductivity_model_lu2007(kay_upper,T(t-1,i),T_std,dug_VWC(t-1,VWC_depth_indices(i)),theta_k,m,Soil_RH, material);
                rho = rho_dry + rho_H2O*dug_VWC(t-1,VWC_depth_indices(i)); %H2O dep
                Cp = tima_specific_heat_model_hillel(rho_dry,rho,dug_VWC(t-1,VWC_depth_indices(i)));
                            %Subsurface Multip Factors (Kieffer, 2013)
                F1 = 2*dt*k(i)/((Cp*rho*layer_size(i)^2)*(1+layer_size(i+1)/layer_size(i)*k(i)/k(i+1)));
                F3 = (1+(layer_size(i+1)/layer_size(i))*(k(i)/k(i+1)))/(1+(layer_size(i-1)/layer_size(i))*(k(i)/k(i-1)));
                F2 = -(1 + F3);
                %Temperature response
                dT = F1*(T(t-1,i+1)+F2*T(t-1,i)+F3*T(t-1,i-1))+dt/(Cp*rho*layer_size(i))*q_evap_i;
                T(t,i) = T(t-1,i)+dT;
            else
                k(i) = kay_lower;
                rho = rho_lower;%1500; %kg/m^3
                Cp = Cp_lower;%2000; %J/kgK
                q_evap_i = 0; %Evap_Coeff
            end
        end
    end
    T_Surf_C = T(:,1) - 273.15;

else %Mapping mode
    % ***************    

    % ***************
    % Loop through observation timesteps
    for t = 2:length(air_temp_C)
        %***************DENSITY******************
        rho = rho_dry + rho_H2O*dug_VWC(t-1,1); %Dry density + water content
        %****************************************
        
        %*********LATENT HEAT***********
        [q_evap_1,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,1),dug_VWC(t-1,1));
        %*******************************

        %*********THERMAL CONDUCTIVITY***********
        k(1) = tima_conductivity_model_lu2007(kay_top,T(t-1,1),T_std,dug_VWC(t-1,1),theta_k,m,Soil_RH, material);%Top
        k(NLAY) = kay_bottom;%Bottom
        %****************************************
        
        %*********Specific heat***********
        Cp = tima_specific_heat_model_hillel(rho_dry,rho,dug_VWC(t-1,1));
        %****************************************

        %*********SENSIBLE HEAT***********
        q_conv = tima_sensible_heat_model(CH,windspeed_horiz(t-1),air_temp_K(t-1),T(t-1,1));
        %*******************************
              
               
        %*********Radiative Flux***********
        phi = solar_azimuth_cwfromS(t-1);
        Z = solar_zenith_apparent(t-1);
        Incidence = cosd(Z)*cosd(slope_angle)+sind(Z)*sind(slope_angle)*cosd(phi-aspect_cwfromS);Incidence(Incidence>1) = 1; Incidence(Incidence<-1) = -1;
        q_rad_full = (1-albedo)*Incidence*(1-f_diff)*r_short_upper(t-1)/cosd(Z)+(1-albedo)*omega*f_diff*r_short_upper(t-1)+(1-albedo)*(1-omega)*r_short_lower(t-1)+omega*emissivity*r_long_upper(t-1)-omega*emissivity*sigma*T(t-1,1)^4; %net heat flux entering surface
        Lit_Fraction = shadow_data(shadow_time_ind(t-1)); %Assumes shadow = 0, and unshadowed = 1;
        q_rad = Lit_Fraction*q_rad_full+(1-Lit_Fraction)*((1-albedo)*omega*f_diff*r_short_upper(t-1)+(1-albedo)*(1-omega)*r_short_lower(t-1)+omega*emissivity*r_long_upper(t-1)-omega*emissivity*sigma*T(t-1,1)^4); %net heat flux entering surface
        %*******************************
        
        %*********Ground Flux***********
        %ref: Kieffer, 2013
        q_G = -k(1)*(T(t-1,1)-T(t-1,2))/layer_size(1); %heat flux from conducting with lower layer
        %*******************************

        %*********Combine Heat transfer elements***********
        W = q_G+q_conv+q_rad+q_evap_1; %Heat Flux, W/m^2
        dT = dt/(Cp*rho*layer_size(1))*(W); %Temperature change for given thickness of material with known volumetric Cp
        T(t,1) = T(t-1,1) + dT; %new T at surface
        %*******************************

        %*********Subsurface Flux (from KRC, Kieffer, 2013)***********
        for i = 2:NLAY-1 %layer loop
%           k_H2O = 1.1e-05*(T(t-1,1)- 273.15).^2 + 0.00234*(T(t-1,1)- 273.15) + 0.552; %Bristow 2002, temp dep
            if i < D_i
                if sum(layer_size(1:i)) <= 0.05 %From Tran 2016, van de Griend and Owe 1994, Camillo and Gurney 1986.
                    [q_evap_i,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,i),dug_VWC(t-1,VWC_depth_indices(i)));
                else
                    [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,i),dug_VWC(t-1,VWC_depth_indices(i)));
                    q_evap_i = 0; %Evap_Coeff
                end
                k(i) = tima_conductivity_model_lu2007(kay_upper,T(t-1,i),T_std,dug_VWC(t-1,VWC_depth_indices(i)),theta_k,m,Soil_RH, material);
                rho = rho_dry + rho_H2O*dug_VWC(t-1,VWC_depth_indices(i)); %H2O dep
                Cp = tima_specific_heat_model_hillel(rho_dry,rho,dug_VWC(t-1,VWC_depth_indices(i)));
            else
                if material == 'ice'
                    k(i) = kay_lower;
                    rho = rho_lower;%1500; %kg/m^3
                    Cp = Cp_lower;%2000; %J/kgK
                elseif material == 'amorphous'
                else
                    k(i) = tima_conductivity_model_lu2007(kay_lower,T(t-1,i),T_std,dug_VWC(t-1,VWC_depth_indices(i)),theta_k,m,Soil_RH, material);
                end
                q_evap_i = 0; %Evap_Coeff
            end
            %Subsurface Multip Factors (Kieffer, 2013)
            F1 = 2*dt*k(i)/((Cp*rho*layer_size(i)^2)*(1+layer_size(i+1)/layer_size(i)*k(i)/k(i+1)));
            F3 = (1+(layer_size(i+1)/layer_size(i))*(k(i)/k(i+1)))/(1+(layer_size(i-1)/layer_size(i))*(k(i)/k(i-1)));
            F2 = -(1 + F3);
            %Temperature response
            dT = F1*(T(t-1,i+1)+F2*T(t-1,i)+F3*T(t-1,i-1))+dt/(Cp*rho*layer_size(i))*q_evap_i;
            T(t,i) = T(t-1,i)+dT;
        end
        %*******************************
        
    end
    T_Surf_C = T(:,1) - 273.15;
end
end
