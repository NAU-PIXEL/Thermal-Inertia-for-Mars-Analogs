# Thermal Inertia for Mars Analogs (TIMA) Model
 TIMA is a finite element numerical model of soil heat transfer that that considers time-variability in thermal conductivity as a function of temperature, moisture, and composition. The detailed parametrization of thermal conductivity is meant to enable direct translation between soil thermophysical properties on Earth and Mars. This repository is maintained by [Ari Koeppel](https://earthsciences.dartmouth.edu/people/ari-koeppel). Please cite Koeppel et al., (2024) when using the model.

# References:
Koeppel, A.H., Edwards, C.S., Edgar, L.A., Nowicki, S., Bennett, K.A., Gullikson, A., Piqueux, S., Eifert, H., Chapline, D. and Rogers, A.D., 2024. A novel surface energy balance method for thermal inertia studies of terrestrial analogs. Earth and Space Science, 11(9), p.e2023EA003259.

# Overview: 
The original implementation of the model was to derive a characteristic soil thermal conductivity under dry conditions at 300K. The procedure achieves this by inputing micrometeorological data into a surface energy balance forward model and adjusting thermal conductivity (along with 5 other modifying parameters) to fit surface skin temperature data (typically obtained through radiometer observation).

# Forward Model:
  [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    windspeed_horiz,T_deep,initial_temps,layer_size,dug_VWC,evap_depth,RH,emissivity,...
    pressure_air_pa,'slope_angle',slope_angle,'aspect_cwfromS','solar_azimuth_cwfromS',solar_azimuth_cwfromS,...
    'solar_zenith_apparent',solar_zenith_apparent,'f_diff',f_diff,'shadow_data',shadow_data,...
    'shadow_time_ind',shadow_time_ind,'MappingMode',true)

# Fitting Process Input/Output Setup:
  TData: Time data - struct of timeseries data variables 
      TData.air_Temp_C: 
      TData.DF:
      TData.dug_VWC_smooth: 
      TData.dug_VWC_smooth_II: 
      TData.evap_depth 
      TData.evap_depth_II 
      TData.humidity: as a fraction
      TData.pressure_air_Pa: 
      TData.r_long_upper:
      TData.r_short_upper: 
      TData.r_short_lower: 
      TData.solarazimuth_cwfromS:
      TData.solarzenith_apparent:
      TData.timed_albedo: 
      TData.timed_albedo_II: 
      TData.TIMESTAMP: 
      TData.temps_to_fit: 
      TData.temps_to_fit_II: 
      TData.windspeed_horiz_ms: 
  MData: Model Data - Struct of static and model format variables
      MData.burnin: 
      MData.dt:
      MData.density: 
      MData.emissivity: 
      MData.layer_size: 
      MData.mccount:  NwalkersxNsteps This is the total number of MC proposals, -NOT the number per chain.
      MData.minit:
      MData.nvars: 
      Mdata.parallel: t or f   [f]
      MData.T_deep: 
      MData.T_start: 
      MData.T_std: 
      MData.ThinChain: [20]
      MData.VWC_depth_indices: 
      MData.notes:
      MData.vars_init:
  outDIR:

# Run Fitting Process
[models,names] = tima_TI_Earth_Tower(TData,MData,outDIR)

# Visualize
tima_plot_results(TData,MData,models,names,varargin)

## Input data:
Example input files are located in the [`in`](hhttps://github.com/NAU-PIXEL/Thermal-Inertia-for-Mars-Analogs/tree/main/Example%20Data) folder.
