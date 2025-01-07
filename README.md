# Thermal Inertia for Mars Analogs (TIMA) Model
 TIMA is a finite element numerical model of soil heat transfer that that considers time-variability in thermal conductivity as a function of temperature, moisture, and composition. The detailed parametrization of thermal conductivity is meant to enable direct translation between soil thermophysical properties on Earth and Mars. This repository is maintained by [Ari Koeppel](https://earthsciences.dartmouth.edu/people/ari-koeppel). Please cite Koeppel et al., (2024) when using the model.

# References:
Koeppel, A.H., Edwards, C.S., Edgar, L.A., Nowicki, S., Bennett, K.A., Gullikson, A., Piqueux, S., Eifert, H., Chapline, D. and Rogers, A.D., 2024. A novel surface energy balance method for thermal inertia studies of terrestrial analogs. Earth and Space Science, 11(9), p.e2023EA003259.

# Overview: 
The original implementation of the model was designed to derive a characteristic dry soil thermal conductivity at 300K. The procedure achieves this by inputing micrometeorological data into a surface energy balance forward model and adjusting thermal conductivity (along with 5 other modifying parameters) to fit surface skin temperature data (typically obtained through radiometer observation). Model updates have included multilayer paramtrizations, including consideration of subsurface ice. Thus, the model can be used to derive soil physical properties, layering and the depth to subsurface material transitions based entirely on surface temperature and weather observations.

# Forward Model:
  [T_surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,
    rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,
    windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,
    pressure_air_pa,varargin)

# Forward Model Inputs:
## Input Parameters ('single quotes' means optional varargin)
  ### Timeseries data
      air_Temp_C: [C] near surface air temperature, typically at 3 m AGL.
          (1D vector)
      'albedo': [0-1] time variant albedo (e.g., due to wetting) for
          fitting surface. (1D vector)
      evap_depth: [m] Depth of the evaporation front (e.g., as function
          of VWC). (1D vector)
      'f_diff': [0-1] Fraction of  Global Horizontal Irradiance (GHI)
          or r_short_upper that is diffuse. (1D vector)
      RH: [0-1] array of near surface relative humidity values,
          typically at 3m AGL. (1D vector)
      pressure_air_pa: [Pa] station pressure, typically at 3m AGL. (1D vector)
      r_long_upper: [W/m^2] Integrated longwave radiation  (4.5 to 42 μm) incident on flat
          surface.(1D vector)
      r_short_lower: [W/m^2] Integrated upwelling shortwave radiation (305 to 2800 nm) from flat
          surface. (1D vector)
      r_short_upper: [W/m^2] Integrated shortwave radiation (305 to 2800 nm) incident on flat
          surface. (1D vector)
      'shadow_time_ind': Index of shadow_data corresponding to
          timeseries. (1D vector)
      'solar_azimuth_cwfromS': [degrees] Solar azimuth in degrees
          clockwise from South, typically -180:180 or 0:360 (1D vector)
      'solar_zenith_apparent': [degrees] Solar zenith in degrees,
          corrected for atmospheric refraction. (1D vector)
      VWC_column: [0-1, decimal fraction by volume] array of volumetric water content
          for each model layer with each time step, typically
          interpolated. (2D vector)
      windspeed_horiz: [m/s] Near surface horizontal wind speed,
          typically at 3m AGL. (1D vector)

  ### Constants:
      'aspect_cwfromS': [degrees] slope aspect clockwise from South. (scalar)
      CE: [Unitless] resistance to latent heat flux coefficient, similar
          to the aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2) (scalar)
      CH: [Unitless] resistance to sensible heat flux coefficient,
          similar to the aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2)
          1/CH should be 0.0028-0.0075 (or CH~100-400) for smooth to
          roughly open soils on Davenport Scale, CH is larger (more
          resistence) with rougher topography or larger vegetation (scalar)
      'depth_transition': [m] Depth of major transition between upper
          material and lower material. (scalar)
      dt: [s] Time step (scalar)
      emissivity: [0-1] Weighted thermal emissivity over wavelength
          range of sensor. (scalar)
      'e_fxn': spectral emissivity function that determines emissivity as
          function of wavelength in um. (function)
      initial_temps: [K] Initial condition, list of center temperatures
          for each layer at start of simulation (1D vector)
      'k_dry_std_lower': [W/mK] bulk dry thermal conductivty of lower layer at T_std (scalar)
      'k_dry_std_mantle': [W/mK] bulk dry thermal conductivty of topmost
          mantling layer at T_std (scalar)
      k_dry_std_upper: [W/mK] bulk dry thermal conductivty of upper layer at T_std (scalar)
      layer_size: [m] List of vertical thickness for each layer
          from top to bottom. (1D vector)
      m: [unitless] pore-network-connectivity-parameter, typically 0-1.3 (scalar)
      'mantle_thickness': [m] thickness of topmost mantling layer. (scalar)
      'material': ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
      'material_lower':  ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at depth (char)
      NDAYS: [unitless] # of days to run equilib model for: More = closer to equilib, fewer = faster (scalar)
      'rho_dry_lower': [kg/m^3] Value for density of lower layer soil.  (scalar)
      rho_dry_upper: [kg/m^3] Value for density of soil.  (scalar)
      'shadow_data': [0-1] Fraction of ROI/pixel in shadow corresponding to time. (1D Vector)
      'slope_angle': [degrees] angle of slope. (scalar)
      theta_k: [0-1, fraction by volume] conductivity soil moisture inflection
          point - in theory this should be similar to saturation value or porosity (scalar) 
      theta_E: [0-1, fraction by volume] latent heat soil moisture inflection
          point (scalar)
      T_deep: [K] Lower boundary condition, fixed temperature (scalar)
      T_std: [K] Standard temperature; typically 300 (scalar)
      'T_adj1': [index, temperature K] pair used to force column
          temperature change at a given time point due to wetting
          (1D vector)
      'T_adj2': [index, temperature K] pair used to force a second
          column temperature change at a given time point due to wetting
          (1D vector)

# Forward Model Output Parameters:
  T_Surf_C = [C] Surface temperature time series (1D vector)

# Fitting Process Input/Output Setup:
  ## Time data - struct of timeseries data variables
      TData.air_Temp_C=Air_Temp_C;
      TData.DF=f_diff;
      TData.VWC_column=VWC_column;
      TData.evap_depth=evap_depth.*ones(size(Air_Temp_C));
      TData.humidity=Humidity;
      TData.humidity=Humidity;
      TData.pressure_air_Pa=Pressure_air_Pa;
      TData.r_long_upper=R_Long_Upper;
      TData.r_short_upper=R_Short_Upper;
      TData.r_short_lower=R_Short_Lower;
      TData.solarazimuth_cwfromS=SolarAzimuthCwfromS;
      TData.solarzenith_apparent=SolarZenith_Apparent;
      TData.timed_albedo=Albedo;
      TData.TIMESTAMP = Data.TIMESTAMP;
      TData.temps_to_fit=Temps_to_fit;
      TData.windspeed_horiz_ms=WindSpeed_ms_10;
      TData.temp_column = Dug_Temp;

   ## Model Data - Struct of static and model format variables
      MData.burnin_fit=burnin_fit;
      MData.burnin_mcmc=burnin_mcmc;
      MData.dt=dt;  
      MData.density=density; 
      MData.emissivity=emissivity;
      MData.fit_ind = fit_ind;
      MData.layer_size=Layer_size_B;
      MData.material=material;
      MData.material_lower=material_lower;
      MData.nwalkers = nwalkers;
      MData.nstep = nstep;
      MData.ndays = NDAYS;
      MData.minit = minit;
      MData.nvars=nvars;
      MData.parallel=true;
      MData.T_deep= T_Deep; 
      MData.T_start= T_Start;
      MData.T_std=T_std;
      MData.ThinChain=20; %[20]
      MData.notes = input('Please enter label/note: ', 's');
      MData.vars_init = Vars_init;
      MData.erf = erf;
      MData.mantle_thickness = mantle_thickness;
      MData.k_dry_std_mantle = k_dry_std_mantle;
      
  outDIR: 'output directory path'

# Run Fitting Process
[models,names] = tima_TI_Earth_Tower(TData,MData,outDIR)

# Visualize
tima_plot_results(TData,MData,models,names,varargin)

## Input data:
Example input files are located in the [`Example Data`](hhttps://github.com/NAU-PIXEL/Thermal-Inertia-for-Mars-Analogs/tree/main/Example%20Data) folder.
