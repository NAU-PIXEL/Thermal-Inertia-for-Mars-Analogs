function [] = tima_plot_results(TData,MData,models,names,varargin)
%% TIMA_PLOT_RESULTS
%   Script to plot MCMC histograms and visualize energy fluxes.
%
% Syntax
%   tima_plot_results(TData,MData,models,names,varargin)
%
% Description
%   Script uses the measured surface temperature for 1 day on repeat as the forcing and then distrubutes heat
%   This function uses observational data and assigned thermophysical properties 
%   to estimate the the temperature profile of the subsurface at the start of the simulation.
%   Assumes homogenous bulk density and bulk dry thermal conductivity with depth
%
% varargin:
%   'TwoSpot': Option to plot results for two different sets of surface temperatures
%       with the same parameters (default=false) (logical) 
%
% Input Parameters
%   TData: Time data - struct of timeseries data variables (all vectors)
%       TData.air_Temp_C: [C] near surface air temperature, typically at 3 m AGL
%       TData.DF: [decimal fraction] Fraction of  Global Horizontal Irradiance (GHI)
%           or r_short_upper that is diffuse (see
%           https://github.com/sandialabs/MATLAB_PV_LIB)
%       TData.temps_to_fit: [C] Surface temperature values to be used for
%           fitting.
%       TData.temp_column: [C] array of Temperatures measured at depth for
%           comparison to model results
%       TData.timed_albedo: [decimal fraction] time variant albedo (e.g.,
%           due to wetting) for fitting surface.
%       TData.TIMESTAMP: datetime array associated with each table row
%           spaced by Mdata.dt [s] (datetime)
%       TData.humidity: [decimal fraction] array of near surface relative humidity values,
%           typically at 3m AGL TData.r_long_upper: [W/m^2] Integrated longwave radiation (4.5 to 42 μm) incident on flat
%           surface
%       TData.pressure_air_Pa: [Pa] station pressure, typically at 3m AGL
%       TData.r_short_lower: [W/m^2] Integrated upwelling shortwave radiation (305 to 2800 nm) from flat
%           surface
%       TData.r_short_upper: [W/m^2] Integrated shortwave radiation (305 to 2800 nm) incident on flat
%           surface
%       TData.solarazimuth_cwfromS: [degrees] Solar azimuth in degrees
%           clockwise from South, typically -180:180 or 0:360
%       TData.solarzenith_apparent: [degrees] Solar zenith in degrees,
%           corrected for atmospheric refraction.
%       TData.VWC_column: [decimal fraction by volume] array of volumetric water content
%           for each model layer with each time step, typically
%           interpolated.
%       TData.windspeed_horiz_ms: [m/s] Near surface horizontal wind speed,
%           typically at 3m AGL.
%   
%   MData: Model Data - Struct of static and model format variables
%       MData.burnin_fit: initial time length (s) to ignore in fitting (vector, default=0)
%       MData.density: [kg/m^3] Value for density of soil beneath tower.  (vector)
%       MData.dt: [s] Time step (vector)
%       MData.emissivity: [0-1] Weighted thermal emissivity over wavelength
%           range of sensor. (vector)
%       MData.erf: Uncertainty as function of observed temperature (function_handle)
%       MData.fit_ind: Indecies of temps_to_fit in which to apply fitting
%           to (vector)
%       MData.layer_size: [m] List of vertical thickness for each layer
%           from top to bottom. (vector)
%       MData.material: ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
%       MData.material_lower:  ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at depth (char)
%       MData.nstep: Number of iterations for curve fitting, 250 is good
%           (vector)
%       MData.nvars: Number of variables being fit for (vector)
%       MData.nwalkers: Number of walkers in MCMC ensemble (vector)
%       MData.T_adj1: [index, temperature K] pair used to force column
%           temperature change at a given time point due to wetting (vector,
%           optional)
%       MData.T_adj2: [index, temperature K] pair used to force a second
%           column temperature change at a given time point due to wetting (vector,
%           optional)
%       MData.T_deep: [K] Lower boundary condition, fixed temperature (vector)
%       MData.T_start: [K] Initial condition, list of center temperatures
%           for each layer at start of simulation (vector)
%       MData.T_std: [K] Standard temperature; typically 300 (vector)
%       MData.vars_init: [k-upper [W/mK], Pore network con. par. (mk) [unitless],...
%           Surf. ex. coef. (CH) [unitless], Surf. ex. coef. (CE) [unitless], Soil Moist. Infl. (thetak) [% by volume],...
%           Soil Moist. Infl. (thetaE) [% by volume], (Transition Depth [m]), (k-lower [W/mK])]
%           List of 6-8 inputs for variables to serve as either initial or fixed values. (vector)
%
%   models: A nvars by (nwalkers*nsteps/ThinChain) matrix with the thinned markov chains (vector)
%
% Author
%    Ari Koeppel -- Copyright 2024


p = inputParser;
p.addRequired('TData',@isstruct);
p.addRequired('MData',@isstruct);
p.addRequired('models');
p.addRequired('names',@iscellstr);
p.addParameter('TwoSpot',false,@islogical);
p.parse(TData, MData, models, names,varargin{:});
p=p.Results;
TwoSpot = p.TwoSpot;
MData.fit_ind = MData.fit_ind(MData.burnin_fit/MData.dt:end);

%% Make corner Plot
% ***************
% NAMED PARAMETERS:
%   range: Restrict visual limits to central [99.5] percentile.
%   names: A cell of strings with labels for each parameter in m.
%   ks: enable kernel smoothing density instead of histograms [false]
%   support: a 2xM matrix with upper and lower limits.
%   ess: effective sample size. [default: auto-determine ess using EACORR.]
%        - used to adjust bandwidth estimates in kernel density estimates.
%   scatter: show scatter plot instead of 2d-kernel density estimate [true if #points<2000]. 
%   fullmatrix: display upper corner of plotmatrix. [false]
%   color: A color-theme for the plot. [.5 .5 .5].
%   grid: show grid. [false].
% ***************
figure
olive = [110/255,117/255,14/255];
[H,RESULTS] = tima_ecornerplot(models,'color',olive,'names',names,'grid',true,'ks',true);


%% Plot Comparison
if size(RESULTS,2) == 6
    formod = @(FitVar) tima_heat_transfer(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,'albedo',TData.timed_albedo,'material',MData.material);%,'e_fxn',MData.e_fxn);
    
    formod_fluxes = @(FitVar) tima_heat_transfer_energy_terms(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,'albedo',TData.timed_albedo,'material',MData.material);%'e_fxn',MData.e_fxn);
     [Result_Temp_Surf,Result_Temp_Sub,q_latent,Result_keff,q_conv,q_rad,q_G] = formod_fluxes(RESULTS(2,:));
elseif size(RESULTS,2) == 7
     formod = @(FitVar) tima_heat_transfer(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,'albedo',TData.timed_albedo,'material',MData.material,'depth_transition',FitVar(7),'material_lower',MData.material_lower);%,'e_fxn',MData.e_fxn);
     
     formod_fluxes = @(FitVar) tima_heat_transfer_energy_terms(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,'albedo',TData.timed_albedo,'material',MData.material,'depth_transition',FitVar(7),'material_lower',MData.material_lower);%'e_fxn',MData.e_fxn);
     [Result_Temp_Surf,Result_Temp_Sub,q_latent,Result_keff,q_conv,q_rad,q_G] = formod_fluxes(RESULTS(2,:));
elseif size(RESULTS,2) == 8
         formod = @(FitVar) tima_heat_transfer(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,'albedo',TData.timed_albedo,'material',MData.material,'depth_transition',FitVar(7),'material_lower',MData.material_lower,'k_dry_std_lower',FitVar(8));%,'e_fxn',MData.e_fxn);
     
     formod_fluxes = @(FitVar) tima_heat_transfer_energy_terms(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,'albedo',TData.timed_albedo,'material',MData.material,'depth_transition',FitVar(7),'material_lower',MData.material_lower,'k_dry_std_lower',FitVar(8));%'e_fxn',MData.e_fxn);
     [Result_Temp_Surf,Result_Temp_Sub,q_latent,Result_keff,q_conv,q_rad,q_G] = formod_fluxes(RESULTS(2,:));
end

if TwoSpot == true
    formod_fluxes = @(FitVar) tima_heat_transfer_energy_terms(FitVar(1),FitVar(2),FitVar(3),...
            FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
            TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
            TData.VWC_II_column,TData.evap_depth_II,TData.humidity,MData.emissivity,...
            TData.pressure_air_Pa,'T_adj1',MData.T_adj1,'T_adj2',MData.T_adj2,'albedo',TData.timed_albedo_II);
    [Result_Temp_Surf,Result_Temp_Sub,q_latent,Result_keff,q_conv,q_rad,q_G] = formod_fluxes(RESULTS(2,:));
end

figure
hold on
xlabel('Time (hr)');
ylabel('Temperature (C)');
if TwoSpot == true
    F(1) = fill([TData.TIMESTAMP(MData.fit_ind); flipud(TData.TIMESTAMP(MData.fit_ind))],[TData.temps_to_fit_II(MData.fit_ind)-MData.erf(TData.temps_to_fit_II(MData.fit_ind));flipud(TData.temps_to_fit_II(MData.fit_ind)+MData.erf(TData.temps_to_fit_II(MData.fit_ind)))],[128 193 219]./255,'Linestyle','none','DisplayName','FLIR error');
    F(2) = scatter(TData.TIMESTAMP(MData.fit_ind),TData.temps_to_fit_II(MData.fit_ind),1,'k.','DisplayName','FLIR Surface Observations');
else
    F(1) = fill([TData.TIMESTAMP(MData.fit_ind); flipud(TData.TIMESTAMP(MData.fit_ind))],[TData.temps_to_fit(MData.fit_ind)-MData.erf(TData.temps_to_fit(MData.fit_ind));flipud(TData.temps_to_fit(MData.fit_ind)+MData.erf(TData.temps_to_fit(MData.fit_ind)))],[128 193 219]./255,'Linestyle','none','DisplayName','FLIR error');
    F(2) = scatter(TData.TIMESTAMP(MData.fit_ind),TData.temps_to_fit(MData.fit_ind),1,'k.','DisplayName','FLIR Surface Observations');
end
set(F(1), 'edgecolor', 'none');
set(F(1), 'FaceAlpha', 0.5);
M = plot(TData.TIMESTAMP(MData.fit_ind),Result_Temp_Surf((MData.fit_ind),1),'r', 'LineWidth', 2 ,'DisplayName','Surface Modeled');

hold off
legend([F(2) M], 'Interpreter','none')
fval = tima_fval_chi2v(TData.temps_to_fit(MData.fit_ind),tima_formod_subset(RESULTS(2,:),MData.fit_ind,formod),MData.erf(TData.temps_to_fit(MData.fit_ind)),MData.nvars);
Cp_std = tima_specific_heat_model_hillel(MData.density,MData.density);
TI =  sqrt(RESULTS(2,1)*MData.density*Cp_std);
TIp = sqrt(RESULTS(3,1)*MData.density*Cp_std);
TIm = sqrt(RESULTS(1,1)*MData.density*Cp_std);
ttl = sprintf('TI Top [Jm^{-2}K^{-1}s^{-12}] = %0.2f, chi_v = %0.2f',TI,fval);%Calculate TI from results
title(ttl)

residuals = TData.temps_to_fit(MData.fit_ind)-tima_formod_subset(RESULTS(2,:),MData.fit_ind,formod);
figure
plot(TData.TIMESTAMP(MData.fit_ind),residuals)
title('residuals')

figure
plot(TData.TIMESTAMP(MData.fit_ind),movmean(gradient(TData.temps_to_fit(MData.fit_ind)),60),'DisplayName','Observed')
hold on
plot(TData.TIMESTAMP(MData.fit_ind),movmean(gradient(tima_formod_subset(RESULTS(2,:),MData.fit_ind,formod)),60),'DisplayName','Modeled')
title('1hr gradients')
legend

%% Plot sub fluxes!
figure
hold on
ylabel('Temperature (C)');
if TwoSpot == true
    plot(TData.TIMESTAMP(MData.fit_ind),TData.temps_to_fit_II(MData.fit_ind),'k','LineWidth', 0.5,'DisplayName','FLIR Surface Observations');
else
    plot(TData.TIMESTAMP(MData.fit_ind),TData.temps_to_fit(MData.fit_ind),'k','LineWidth', 0.5,'DisplayName','FLIR Surface Observations');
end
title('FLIR Surface Observations')


figure
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP(MData.fit_ind),q_conv(MData.fit_ind),'g', 'LineWidth', 1,'DisplayName','sensible heat');
title('Sensible')

figure
hold on
ylabel('W/mK');
plot(TData.TIMESTAMP(MData.fit_ind),Result_keff((MData.fit_ind),1),'r:', 'LineWidth', 1 ,'DisplayName','k_eff (layer 1)');
plot(TData.TIMESTAMP(MData.fit_ind),Result_keff((MData.fit_ind),2),'r', 'LineWidth', 1 ,'DisplayName','k_eff (layer 2)');
plot(TData.TIMESTAMP(MData.fit_ind),Result_keff((MData.fit_ind),end-1),'r--', 'LineWidth', 1 ,'DisplayName','k_eff (layer end-1)');
hold off
legend
title('K_eff')

figure
hold on
ylabel('W/mK');
plot(TData.TIMESTAMP(MData.fit_ind),Result_Temp_Sub((MData.fit_ind),1),'k:', 'LineWidth', 1 ,'DisplayName','Layer 1');
plot(TData.TIMESTAMP(MData.fit_ind),Result_Temp_Sub((MData.fit_ind),2),'k', 'LineWidth', 1 ,'DisplayName','Layer 2');
plot(TData.TIMESTAMP(MData.fit_ind),Result_Temp_Sub((MData.fit_ind),3:end-3),'b', 'LineWidth', 0.25 ,'DisplayName','Middle Layers');
plot(TData.TIMESTAMP(MData.fit_ind),TData.temp_column((MData.fit_ind),:),'g', 'LineWidth', 0.25 ,'DisplayName','Measured');
plot(TData.TIMESTAMP(MData.fit_ind),Result_Temp_Sub((MData.fit_ind),end-2),'k-.', 'LineWidth', 1 ,'DisplayName','Layer end-2');
plot(TData.TIMESTAMP(MData.fit_ind),Result_Temp_Sub((MData.fit_ind),end-1),'k--', 'LineWidth', 1 ,'DisplayName','Layer end-1');
hold off
legend
title('Temperature')

figure
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP(MData.fit_ind),q_rad(MData.fit_ind),'m', 'LineWidth', 1 ,'DisplayName','Radiative heat');
title('Radiative')


figure
for time = 1:length(TData.temps_to_fit)
    full_latent(time) = sum(q_latent(time,:));
end
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP(MData.fit_ind),full_latent(MData.fit_ind),'c', 'LineWidth', 1 ,'DisplayName','latent heat');
title('Latent')

figure
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP(MData.fit_ind),q_G(MData.fit_ind),'b', 'LineWidth', 1,'DisplayName','ground heat');
title('Ground')
end