function [T_Surf_C] = tima_heat_transfer_skin_treatment(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,evap_depth,RH,emissivity,...
    pressure_air_pa,varargin)
%***************
% TIMA_HEAT_TRANSFER
%   This function uses observational data and assigned thermophysical properties 
%   to estimate the emitting surface temperature through time.

% Syntax
%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,evap_depth,RH,emissivity,...
    % pressure_air_pa)

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,evap_depth,RH,emissivity,...
    % pressure_air_pa,'material',material)

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,evap_depth,RH,emissivity,...
    % pressure_air_pa,'T_adj1',[index, T_K])

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,evap_depth,RH,emissivity,...
    % pressure_air_pa,k_dry_std_lower,depth_transition)

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,evap_depth,RH,emissivity,...
    % pressure_air_pa,k_dry_std_lower,depth_transition,'rho_dry_lower',rho_dry_lower)

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,evap_depth,RH,emissivity,...
    % pressure_air_pa,'albedo',albedo)

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,evap_depth,RH,emissivity,...
    % pressure_air_pa,'slope_angle',slope_angle,'aspect_cwfromS','solar_azimuth_cwfromS',solar_azimuth_cwfromS,...
    % 'solar_zenith_apparent',solar_zenith_apparent,'f_diff',f_diff,'shadow_data',shadow_data,...
    % 'shadow_time_ind',shadow_time_ind,'MappingMode',true)

% Description
%   This model uses a surrogate optimization solver with 1-7 paramters to fit IR-image-observed
%   surface temperatures using meteorological data and surface parameters derived from
%   aerial/satellite imagery. Fitting Parameters can include any or all of: top layer thermal
%   conductivity at 300K, bottome layer thermal conductivity at 300K, Depth of thermal conductivity
%   transition, pore-network-connectivity, Surf. ex. coef. (sensible), Surf. ex. coef. (latent), Soil Moist. Inflection (conductivity) (%), Soil Moist. Infl. (latent heat) (%)

% Input Parameters
%   Mapping Mode: True or False
%   Observational data
%       dt
%       Air_Temp_C
%       R_Short_Upper: must be positive
%       R_Short_Lower: Must be positive
%       R_Long_Upper
%       WindSpeed_horiz
%       Dug_VWC
%       evap_depth
%       VWC_depth_indices
%       RH
%       rho_dry
%       emissivity
%       albedo
%       pressure_air
%       f_diff
%       material
%    Boundary conditions
%       T_std
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

%   OPTIONAL
%       T_adj1: [index, T_K] pair used to adjust temperature during wetting
%       T_adj2: [index, T_K] pair used to adjust temperature during wetting

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
p.addRequired('r_short_lower',@(x)all(x>=0) && isnumeric(x));
p.addRequired('r_short_upper',@(x)all(x>=0) && isnumeric(x));
p.addOptional('k_dry_std_lower',k_dry_std_upper,@isnumeric);
p.addOptional('depth_transition',sum(layer_size),@isnumeric);


p.addParameter('rho_dry_lower',rho_dry_upper,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('albedo',[],@isnumeric);
p.addParameter('slope_angle',0,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('aspect_cwfromS',[],@(x)isnumeric(x)&&isscalar(x));
p.addParameter('solar_azimuth_cwfromS',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('solar_zenith_apparent',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('f_diff',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('shadow_data',[],@isnumeric);
p.addParameter('shadow_time_ind',[],@isnumeric);
p.addParameter('material',"basalt",@ischar);
p.addParameter('MappingMode',false,@islogical);
p.addParameter('T_adj1',[],@(x)isnumeric(x)&&length(x)==2);
p.addParameter('T_adj2',[]);
p.parse(r_short_lower, r_short_upper, varargin{:});

p=p.Results;

k_dry_std_lower = p.k_dry_std_lower;
depth_transition = p.depth_transition;

rho_dry_lower = p.rho_dry_lower;
albedo = p.albedo;
slope_angle = p.slope_angle;
aspect_cwfromS = p.aspect_cwfromS;
solar_azimuth_cwfromS = p.solar_azimuth_cwfromS;
solar_zenith_apparent = p.solar_zenith_apparent;
f_diff = p.f_diff;
material = p.material;
MappingMode = p.MappingMode;
T_adj1 = p.T_adj1;if ~isempty(T_adj1), T_adj1(1) = T_adj1(1)/(dt/60);end
T_adj2 = p.T_adj2;if ~isempty(T_adj2), T_adj2(1) = T_adj2(1)/(dt/60);end

% ***************

% ***************
% Inputs and constants:
air_temp_K = air_temp_C+273.15;
% k_H2O = 0.6096; %from Haigh 2012
rho_H2O = 997; %kg/m^3, does not vary much in T range so assume constant density
sigma = 5.670374419e-8; % Stefan-Boltzmann Const
omega = (1+cosd(slope_angle))/2; % visible sky fraction ISOtropic sky model
% T_std = 300; %Standard temperature K
% ***************
% Initial and boundary conditions setup + array initiation:
NLAY = length(layer_size);
% T = NaN(NLAY,1); %Matrix of Temperatures for each layer and time point, T(1,t) is surface temp array that will be fit to obs data
T_surf = NaN(length(air_temp_C),1);
T_Surf_C = T_surf;
T_surf(1) = initial_temps(1);
T = initial_temps;
T(NLAY) = T_deep; %Constant lower boundary temp
k = k_dry_std_upper.*ones(1,NLAY);
kay_upper = k_dry_std_upper;
kay_lower = k_dry_std_lower;
[~,D_z] = min(abs(depth_transition-cumsum(layer_size))); %index of depth for k transition
FC = NaN(NLAY-1,1);
RLAY = layer_size(2)/layer_size(1);
for z = 1:NLAY-1
    FC(z) = 2*dt/layer_size(z)^2;
end

for t = 1:length(air_temp_C)
    %*********LATENT HEAT & THERMAL CONDUCTIVITY***********
    [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t),windspeed_horiz(t),RH(t),air_temp_K(t),T(NLAY),dug_VWC(t,end));
    k(NLAY) = tima_conductivity_model_lu2007(kay_lower,T(NLAY),T_std,dug_VWC(t,end),theta_k,m,Soil_RH,material);%Top;%Bottom
    [q_evap_1,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t),windspeed_horiz(t),RH(t),air_temp_K(t),T(2),dug_VWC(t,1));
    k(1) = tima_conductivity_model_lu2007(kay_upper,T(2),T_std,dug_VWC(t,1),theta_k,m,Soil_RH,material);%Top
    %****************************************
    
    %*********Specific heat***********
    % rho = rho_dry_upper + rho_H2O*dug_VWC(t,1); %Dry density + water content
    % Cp = tima_specific_heat_model_DV1963(rho_dry_upper,rho,T(1),material);%tima_specific_heat_model_hillel(rho_dry_upper,rho);%
    %****************************************
    GGT = 0.1;
    delta = GGT+1;
    if t>1
        T_surf(t)=T_surf(t-1);
    end
    count = 1;
    while delta > GGT && t > 1 %newton-raphson iteration to solve for ground temperatures    
        %*********SENSIBLE HEAT***********
        [q_conv,SHCoeff] = tima_sensible_heat_model(CH,windspeed_horiz(t),air_temp_K(t),T_surf(t));
        %*******************************   
        %*********Radiative Flux***********
        if MappingMode == false %Tower mode
            if isempty(solar_azimuth_cwfromS) || isempty(solar_zenith_apparent) || isempty(aspect_cwfromS) || isempty(f_diff) && slope_angle == 0
                if isempty(albedo)
                    q_rad = r_short_upper(t)-r_short_lower(t)+emissivity*r_long_upper(t)-emissivity*sigma*T_surf(t)^4;% %net heat flux entering surface assuming no transmission (Yes emissivity term in down)
                else
                    q_rad = (1-albedo(t))*r_short_upper(t)+emissivity*r_long_upper(t)-emissivity*sigma*T_surf(t)^4; %net heat flux entering surface
                end
            else
                phi = solar_azimuth_cwfromS(t);
                Z = solar_zenith_apparent(t);
                Incidence = cosd(Z)*cosd(slope_angle)+sind(Z)*sind(slope_angle)*cosd(phi-aspect_cwfromS);Incidence(Incidence>1) = 1; Incidence(Incidence<-1) = -1;
                q_rad = (1-albedo(t))*Incidence*(1-f_diff)*r_short_upper(t)/cosd(Z)+(1-albedo(t))*omega*f_diff*r_short_upper(t)+(1-albedo(t))*(1-omega)*r_short_lower(t)+omega*emissivity*r_long_upper(t)-omega*emissivity*sigma*T_surf(t)^4; %net heat flux entering surface
            end
        else %Mapping Mode
            phi = solar_azimuth_cwfromS(t);
            Z = solar_zenith_apparent(t);
            Incidence = cosd(Z)*cosd(slope_angle)+sind(Z)*sind(slope_angle)*cosd(phi-aspect_cwfromS);Incidence(Incidence>1) = 1; Incidence(Incidence<-1) = -1;
            q_rad_full = (1-albedo)*Incidence*(1-f_diff)*r_short_upper(t)/cosd(Z)+(1-albedo)*omega*f_diff*r_short_upper(t)+(1-albedo)*(1-omega)*r_short_lower(t)+omega*emissivity*r_long_upper(t)-omega*emissivity*sigma*T_surf(t)^4; %net heat flux entering surface
            Lit_Fraction = shadow_data(shadow_time_ind(t)); %Assumes shadow = 0, and unshadowed = 1;
            q_rad = Lit_Fraction*q_rad_full+(1-Lit_Fraction)*((1-albedo)*omega*f_diff*r_short_upper(t)+(1-albedo)*(1-omega)*r_short_lower(t)+omega*emissivity*r_long_upper(t)-omega*emissivity*sigma*T_surf(t)^4); %net heat flux entering surface
        end
        %*******************************
        
        %*********Ground Flux***********
        %ref: Kieffer, 2013
        q_G_surf = k(1)*(T(2)-T_surf(t))/(layer_size(2)/2); %heat flux from conducting with lower layer
        %*******************************
    
        %*********Combine Heat transfer elements***********
        W = q_G_surf+q_conv+q_rad+q_evap_1; %Heat Flux, W/m^2
        dWdt = 4*emissivity*sigma*T_surf(t)^3+2*k(1)/layer_size(2)+SHCoeff;%+LECoeff*(tima_psaturationpa_Sonntag1994(T_surf(t))-tima_psaturationpa_Sonntag1994(T_surf(t-1)))/dt;
        %leaving off LE because small and mostly not relevent for surface skin unless frost
        delta = W/dWdt;
        % q_rad = q_rad + 4*emissivity*sigma*T_surf(t)^3*delta; %For printing out
        % q_conv = q_conv + SHCoeff*delta; %for printing out
        % q_evap_1 = q_evap_1 + LHCoeff*(tima_psaturationpa_Sonntag1994(T_surf(t))-tima_psaturationpa_Sonntag1994(T_surf(t-1)))/dt*delta;
        % q_G_surf = q_G_surf + 2*k(1)/layer_size(2)*delta; %for printing out
        %*******************************
        T_surf(t) = T_surf(t)+delta; %new T at surface
        %ALTERNATE TAYLOR SERIES;
        %T_surf(t) = 
        if ~isempty(T_adj1)
            if dug_VWC(t,1) > 0.03 && (t == ceil(T_adj1(1)))
                T_surf(t) = T_adj1(2); %Value taken from T109 probe in watering can
                break
            end
        end
        if ~isempty(T_adj2)
            if dug_VWC(t,1) > 0.02 && (t == ceil(T_adj2(1)))
                T_surf(t) = T_adj2(2); %Value taken from T109 probe in watering can
                break
            end
        end
        if abs(delta)/T_surf(t) > 0.8
            sprintf('Non convergance for k: %0.4f, m: %0.4f, CH: %0.4f, CE: %0.4f, theta_k: %0.4f, theta_E: %0.4f',k_dry_std_upper,m,CH,CE,theta_k,theta_E)
            T_Surf_C(:)=666;
            return
        elseif abs(delta)/T_surf(t) > 0.1
            delta = 0.7*delta;
        end
        count = count+1;
        if count == 50
            break
        end
    end
    T(1) = T(2)-(1+1/RLAY)*(T(2)-T_surf(t));
    if ~isempty(T_adj1)
        if dug_VWC(t,1) > 0.03 && (t == ceil(T_adj1(1)))
            T(1) = T_adj1(2); %Value taken from T109 probe in watering can
        end
    end
    if ~isempty(T_adj2)
        if dug_VWC(t,1) > 0.02 && (t == ceil(T_adj2(1)))
            T(1) = T_adj2(2); %Value taken from T109 probe in watering can
        end
    end
    %*********Subsurface Flux (from KRC, Kieffer, 2013)***********
    for z = 2:NLAY-1 %layer loop
        if z < D_z
            if sum(layer_size(1:z)) <= evap_depth(t) %Wang 2016 inflection point
                [q_evap_z,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t),windspeed_horiz(t),RH(t),air_temp_K(t),T(z),dug_VWC(t,z-1));
            else
                [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t),windspeed_horiz(t),RH(t),air_temp_K(t),T(z),dug_VWC(t,z-1));
                q_evap_z = 0; %Evap_Coeff
            end
            k(z) = tima_conductivity_model_lu2007(kay_upper,T(z),T_std,dug_VWC(t,z-1),theta_k,m,Soil_RH,material);
            rho = rho_dry_upper + rho_H2O*dug_VWC(t,z-1); %H2O dep
            Cp = tima_specific_heat_model_DV1963(rho_dry_upper,rho,T(z),material);%tima_specific_heat_model_hillel(rho_dry_upper,rho);%
                        %Subsurface Multip Factors (Kieffer, 2013)
        else
            if material == "ice"
                k(z) = kay_lower;
                rho = 1500; %kg/m^3
                Cp = Cp_lower;%2000; %J/kgK
                T(z) = T_Deep;
                q_evap_z = 0; %Evap_Coeff
                Soil_RH = 1;
            else
                if sum(layer_size(1:z)) <= evap_depth(t)
                    [q_evap_z,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t),windspeed_horiz(t),RH(t),air_temp_K(t),T(z),dug_VWC(t,z-1));
                else
                    [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t),windspeed_horiz(t),RH(t),air_temp_K(t),T(z),dug_VWC(t,z-1));
                    q_evap_z = 0; %Evap_Coeff
                end
                k(z) = tima_conductivity_model_lu2007(kay_lower,T(z),T_std,dug_VWC(t,z-1),theta_k,m,Soil_RH,material);
                rho = rho_dry_lower + rho_H2O*dug_VWC(t,z-1); %H2O dep
                Cp = tima_specific_heat_model_DV1963(rho_dry_upper,rho,T(z),material);%tima_specific_heat_model_hillel(rho_dry_upper,rho);%
            end
        end
        %Subsurface Multip Factors (Kieffer, 2013)
        % if layer_size(z)/sqrt(2*dt*(k(z)/rho/Cp)) < 1
        %     error('Non convergance for k: %0.4f, m: %0.4f, CH: %0.4f, CE: %0.4f, theta_k: %0.4f, theta_E: %0.4f',k_dry_std_upper,m,CH,CE,theta_k,theta_E)
        % end
        F1 = FC(z)*k(z)/((Cp*rho)*(1+RLAY*k(z)/k(z+1)));
        F3 = (1+RLAY*(k(z)/k(z+1)))/(1+(1/RLAY)*(k(z)/k(z-1)));
        F2 = -(1 + F3);
        %Temperature response
        dT = F1*(T(z+1)+F2*T(z)+F3*T(z-1))+dt/(Cp*rho*layer_size(z))*q_evap_z;
        T(z) = T(z)+dT;
        if ~isempty(T_adj1)
            if dug_VWC(t,z-1) > 0.02 && (t == ceil(T_adj1(1)))
                T(z) = T_adj1(2); %Value taken from T109 probe in watering can
            end
        end
        if ~isempty(T_adj2)
            if dug_VWC(t,z-1) > 0.02 && (t == ceil(T_adj2(1)))
                T(z) = T_adj2(2); %Value taken from T109 probe in watering can
            end
        end
    end
end
T_Surf_C = T_surf - 273.15;

