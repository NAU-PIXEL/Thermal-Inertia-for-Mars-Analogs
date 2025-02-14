function [models,names] = tima_TI_Earth_Tower(TData,MData,outDIR,varargin)
%% TIMA_TI_EARTH_TOWER (Thermal Inertia for Mars Analogs)
%   Surface energy balance model for deriving thermal inertia in terrestrial sediments using diurnal
%   observations taken in the field to fit 1D multi-parameter model.
%
% Description
%   This model uses an MCMC probabilistic solver with 6-8 paramters to fit IR-image-observed
%   surface temperatures using meteorological data and surface parameters derived from
%   aerial/satellite imagery. Fitting Parameters can include any or all of: top layer bulk dry thermal
%   conductivity at a standard temperature, pore-network-connectivity parameter, Surf. ex. coef. (sensible),
%   Surf. ex. coef. (latent), Soil Moist. Inflection (conductivity) (%),
%   Soil Moist. Infl. (latent heat) (%), depth transition (m), bottom layer bulk dry thermal
%   conductivity at a standard temperature
%
% Syntax
%   [models,names] = tima_TI_Earth_Tower(TData,MData,outDIR)
%   [models,names] = tima_TI_Earth_Tower(TData,MData,outDIR,varargin)
%   [models,names] = tima_TI_Earth_Tower(TData,MData,outDIR,'Mode','1layer','Parallel',true,'SaveModels',true)
%
% varargin options:
%   'IncreaseSampling': Option to increase sampling rate if model does not
%       converge. Increases compute time, but avoids nonconvergence. (default=false) (logical)
%   'Initialize': Option to reInitialize subsurface temperatures for each run. Increases compute time, 
%       but avoids dramatic shifts with parameter optimization leading to instability. (default=false) (logical)
%   'Mode': Fitting mode (options: '1layer', '2layer', '2layer_fixed_lower','2layer_fixed_depth' 
%       -- 2Layer includes fitting of bulk dry thermal conductivity of lower 
%       layer and the transition depth) (default='1layer') (string)
%   'Parallel': Run in ensemble of walkers in parallel. (default=false)
%       (logical)
%   'SaveModels': Option to save the models output to a .mat file
%       (default=false) (logical)
%
% Input Parameters
%   TData: Time data - struct of timeseries data variables (all vectors)
%       TData.air_Temp_C: [C] near surface air temperature, typically at 3
%           m AGL (1D vector)
%       TData.DF: [decimal fraction] Fraction of  Global Horizontal Irradiance (GHI)
%           or r_short_upper that is diffuse (1D vector)
%       TData.evap_depth: [m] Depth of the evaporation front. (1D vector)
%       TData.temps_to_fit_interp: [C] Surface temperature values to be used for
%           fitting with any data gaps filled in. (1D vector)
%       TData.timed_albedo: [decimal fraction] time variant albedo (e.g.,
%           due to wetting) for fitting surface. (1D vector)
%       TData.TIMESTAMP: datetime array associated with each table row
%           spaced by Mdata.dt [s] (1D datetime)
%       TData.humidity: [decimal fraction] array of near surface relative humidity values,
%           typically at 3m AGL TData.r_long_upper: [W/m^2] Integrated longwave radiation (4.5 to 42 μm) incident on flat
%           surface (1D vector)
%       TData.pressure_air_Pa: [Pa] station pressure, typically at 3m AGL (1D vector)
%       TData.r_long_upper: [W/m^2] Integrated longwave radiation (4500 to 42000 nm) incident on flat
%           surface (1D vector)
%       TData.r_short_lower: [W/m^2] Integrated upwelling shortwave radiation (305 to 2800 nm) from flat
%           surface (1D vector)
%       TData.r_short_upper: [W/m^2] Integrated shortwave radiation (305 to 2800 nm) incident on flat
%           surface (1D vector)
%       TData.solarazimuth_cwfromS: [degrees] Solar azimuth in degrees
%           clockwise from South, typically -180:180 or 0:360 (1D vector)
%       TData.solarzenith_apparent: [degrees] Solar zenith in degrees,
%           corrected for atmospheric refraction. (1D vector)
%       TData.VWC_column: [decimal fraction by volume] array of volumetric water content
%           for each model layer with each time step, typically
%           interpolated. (2D vector)
%       TData.windspeed_horiz_ms: [m/s] Near surface horizontal wind speed,
%           typically at 3m AGL. (1D vector)
%   
%   MData: Model Data - Struct of static and model format variables
%       MData.burnin_fit: initial time length (s) to ignore in fitting (vector, default=0)
%       MData.burnin_mcmc: fraction of the chain that should be removed.
%           (scalar, default=0)
%       MData.density: [kg/m^3] Value for density of soil beneath tower.
%           (scalar)
%       MData.dt: [s] Time step (scalar)
%       MData.emissivity: [0-1] Weighted thermal emissivity over wavelength
%           range of sensor. (scalar)
%       MData.erf: Uncertainty as function of observed temperature (function_handle)
%       MData.fit_ind: Indecies of temps_to_fit_interp in which to apply fitting
%           to (scalar)
%       MData.k_dry_std_mantle: [W/mK] bulk dry thermal conductivty of topmost
%           mantling layer at T_std (scalar)
%       MData.layer_size: [m] List of vertical thickness for each layer
%           from top to bottom. (1D vector)
%       MData.mantle_thickness: [m] thickness of topmost mantling layer. (scalar)
%       MData.material: ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
%       MData.material_lower:  ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at depth (char)
%       MData.minit: Vector of Initialized test variables with as many
%           randomized samples as desired for fitting, 50 is good such that
%           vector is nvarx50 (2D vector)
%       MData.notes: Details to record in data structure (string)
%       MData.nstep: Number of iterations for curve fitting, 250 is good
%           (scalar)
%       MData.nvars: Number of variables being fit for (scalar)
%       MData.nwalkers: Number of walkers in MCMC ensemble (scalar)
%       MData.ThinChain: MCMC data reduction by thinning all output chains by only storing every N'th step (scalar, default=10)
%       MData.T_adj1: [index, temperature K] pair used to force column
%           temperature change at a given time point due to wetting (1D vector,
%           optional)
%       MData.T_adj2: [index, temperature K] pair used to force a second
%           column temperature change at a given time point due to wetting (1D vector,
%           optional)
%       MData.T_deep: [K] Lower boundary condition, fixed temperature (scalar)
%       MData.T_start: [K] Initial condition, list of center temperatures
%           for each layer at start of simulation (scalar)
%       MData.T_std: [K] Standard temperature; typically 300 (scalar)
%       MData.vars_init: [k-upper [W/mK], Pore network con. par. (mk) [unitless],...
%           Surf. ex. coef. (CH) [unitless], Surf. ex. coef. (CE) [unitless], Soil Moist. Infl. (thetak) [% by volume],...
%           Soil Moist. Infl. (thetaE) [% by volume], (Transition Depth [m]), (k-lower [W/mK])]
%           List of 6-8 inputs for variables to serve as either initial or fixed values. (1D vector)
%
%   out_DIR: Full path to directory for outputs (string)
%
% Output Parameters
%   models: A nvars by (nwalkers*nsteps/ThinChain) matrix with the thinned markov chains (2D vector)
%   names: array of names of fitting variables (character array)
%
% Author
%    Ari Koeppel -- Copyright 2023
%
% Sources
%   Surrogate Optimization: https://www.mathworks.com/help/gads/table-for-choosing-a-solver.html
%   MCMC: Goodman & Weare (2010), Ensemble Samplers With Affine Invariance, Comm. App. Math. Comp. Sci., Vol. 5, No. 1, 65–80
%       Foreman-Mackey, Hogg, Lang, Goodman (2013), emcee: The MCMC Hammer, arXiv:1202.3665
%       https://github.com/grinsted/gwmcmc Aslak Grinsted 2015
%   

format shortG
p = inputParser;
p.addRequired('TData',@isstruct);
p.addRequired('MData',@isstruct);
p.addRequired('outDIR',@ischar);
p.addParameter('Mode','1layer',@ischar);
p.addParameter('Parallel',false,@islogical);
p.addParameter('SaveModels',false,@islogical);
p.addParameter('IncreaseSampling',false,@islogical);
p.addParameter('Initialize',false,@islogical);
p.parse(TData, MData, outDIR, varargin{:});
p=p.Results;
Mode = p.Mode;
Parallel= p.Parallel;
SaveModels = p.SaveModels;
MData.fit_ind = MData.fit_ind(ceil(MData.burnin_fit/MData.dt):end);
Initialize = p.Initialize;
IncreaseSampling = p.IncreaseSampling;

if strcmp(Mode,'1layer')
    formod = @(FitVar) tima_full_model(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,TData.temps_to_fit_interp,MData.ndays,...
        'material',MData.material,'mantle_thickness',MData.mantle_thickness,...
        'k_dry_std_mantle',MData.k_dry_std_mantle,'IncreaseSampling',IncreaseSampling,'Initialize',Initialize);
    MData.minit = MData.minit(1:6,:);
    MData.lbound = [0.024 0.01 1 1 0.05 0.01];
    MData.ubound = [3.7 1.3 1000 10000 0.9 0.75];
elseif strcmp(Mode,'2layer')
    formod = @(FitVar) tima_full_model(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,TData.temps_to_fit_interp,MData.ndays,'material',MData.material,'depth_transition',FitVar(7),...
        'k_dry_std_lower',FitVar(8),'material_lower',MData.material_lower,...
            'mantle_thickness',MData.mantle_thickness,'k_dry_std_mantle',MData.k_dry_std_mantle,'IncreaseSampling',IncreaseSampling,'Initialize',Initialize);
    MData.minit = MData.minit(1:8,:);
    MData.lbound = [0.024 0.01 1 1 0.05 0.01 0.024 0];
    MData.ubound = [3.7 1.3 1000 10000 0.9 0.75 3.7 5];
elseif strcmp(Mode,'2layer_fixed_depth') % This won't work with MCMC in the set up in ln prior right now
    formod = @(FitVar) tima_full_model(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,TData.temps_to_fit_interp,MData.ndays,'material',MData.material,'depth_transition',...
        MData.vars_init(7),'k_dry_std_lower',FitVar(7),'material_lower',MData.material_lower,...
            'mantle_thickness',MData.mantle_thickness,'k_dry_std_mantle',MData.k_dry_std_mantle,'IncreaseSampling',IncreaseSampling,'Initialize',Initialize);
    MData.minit = MData.minit([1:6 8],:);
    MData.lbound = [0.024 0.01 1 1 0.05 0.01 0.024];
    MData.ubound = [3.7 1.3 1000 10000 0.9 0.75 3.7];
elseif strcmp(Mode,'2layer_fixed_lower')
    formod = @(FitVar) tima_full_model(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,TData.temps_to_fit_interp,MData.ndays,'material',MData.material,...
        'depth_transition',FitVar(7),'k_dry_std_lower',MData.vars_init(8),'material_lower',MData.material_lower,...
            'mantle_thickness',MData.mantle_thickness,'k_dry_std_mantle',MData.k_dry_std_mantle,'IncreaseSampling',IncreaseSampling,'Initialize',Initialize);
    MData.minit = MData.minit(1:7,:);
    MData.lbound = [0.024 0.01 1 1 0.05 0.01 0];
    MData.ubound = [3.7 1.3 1000 10000 0.9 0.75 5];
else
    error('Mode entered does not match available options.')
end
%% Optimiztion
optimize_MCMC = 1; %Option to use MCMC or Surrogate Opt
if optimize_MCMC == 1
    k_air = 0.024; % Tsilingiris 2008
    k_H2O = 0.61;
    if strcmp(MData.material,'basalt') %amorphous', 'granite', 'sandstone', 'clay', 'salt', 'ice'
        k_solid = 2.2; % Bristow, 2002 - basaslt = 2.2, clay = 2.9, qtz = 8.8, Ice = 2.18, Granite = 2.0.
    elseif strcmp(MData.material,'granite')
        k_solid = 2.0; 
    elseif strcmp(MData.material,'clay')
        k_solid = 2.9; 
    elseif strcmp(MData.material,'sandstone')
        k_solid = 8.8; 
    elseif strcmp(MData.material,'salt')
        k_solid = 2.0; 
    elseif strcmp(MData.material,'amorphous')
        k_solid = 2.0; 
    elseif strcmp(MData.material,'ice')
        k_solid = 2.18; 
    else
        k_solid = 2.2;
    end
    logPfuns = {@(theta)tima_ln_prior(theta,'depth_max',sum(MData.layer_size)) @(theta)tima_logPfun_chi2v(TData.temps_to_fit_interp(MData.fit_ind),tima_formod_subset(theta,MData.fit_ind,formod),MData.erf(TData.temps_to_fit_interp(MData.fit_ind)))};
    
    %% Run emcee
    % ***************
    % Named Parameter-Value pairs:
    %   'StepSize': unit-less stepsize (default=2.5).
    %   'ThinChain': Thin all the chains by only storing every N'th step (default=10)
    %   'ProgressBar': Show a text progress bar (default=true)
    %   'Parallel': Run in ensemble of walkers in parallel. (default=false)
    %   'BurnIn': fraction of the chain that should be removed. (default=0)
    % ***************
    % – rejection rate is how often the new parameters are accepted. If this is far from ~30% (meaning inefficient), change the proposal step size (sigma),
    %e.g., rej rate too big , make step size smaller, if too small, make step size bigger
    tic
    [models,LogPs]=gwmcmc(MData.minit,logPfuns,MData.nwalkers*MData.nstep,'BurnIn',MData.burnin_mcmc,'Parallel',Parallel,'ThinChain',MData.ThinChain);
    elapsedTime = toc
    if (size(models,1)<size(models,2))&&(ismatrix(models)), models=models'; end
    if ndims(models)==3
        models=models(:,:)'; 
    end
    M=size(models,2);
    for r=1:M
        quant=quantile(models(:,r),[0.16,0.5,0.84]);
        RESULTS(:,r) = quant;
    end
    
    fval = tima_fval_chi2v(TData.temps_to_fit_interp(MData.fit_ind),tima_formod_subset(RESULTS(2,:),MData.fit_ind,formod),MData.erf(TData.temps_to_fit_interp(MData.fit_ind)),MData.nvars);
    Cp_std = tima_specific_heat_model_DV1963(MData.density,MData.density,300,MData.material);
    TI =  sqrt(RESULTS(2,1)*MData.density*Cp_std);
    TIp = sqrt(RESULTS(3,1)*MData.density*Cp_std);
    TIm = sqrt(RESULTS(1,1)*MData.density*Cp_std);
else
    opts = optimoptions('surrogateopt','InitialPoints',MData.minit,'UseParallel',Parallel,'MaxFunctionEvaluations',MData.nstep);
    Obj = @(theta) tima_fval_chi2v(TData.temps_to_fit_interp(MData.fit_ind),tima_formod_subset(theta,MData.fit_ind,formod),MData.erf(TData.temps_to_fit_interp(MData.fit_ind)),MData.nvars);
    problem = struct('solver','surrogateopt','lb',MData.lbound,'ub',MData.ubound,'objective',Obj,'options',opts) ; 
    [RESULTS_holder,fval] = surrogateopt(problem);
    RESULTS(:) = RESULTS_holder';
    Cp_std = tima_specific_heat_model_DV1963(MData.density,MData.density,300,MData.material);
    TI =  sqrt(RESULTS(2,1)*MData.density*Cp_std);
    TIp = sqrt(RESULTS(3,1)*MData.density*Cp_std);
    TIm = sqrt(RESULTS(1,1)*MData.density*Cp_std);
end
%% Print Log File 
c = fix(clock);
fname = sprintf('tima_output_%02.0f%02.0f-%s.txt',c(4),c(5),date);
names = {'k-dry-300-upper (W/mK)' 'Pore network con. par. (mk)' 'Surf. ex. coef. (CH)' 'Surf. ex. coef. (CE)' 'Soil Moist. Infl. (thetak) (%)' 'Soil Moist. Infl. (thetaE) (%)' 'Depth Tansition (m)' 'k-dry-300-lower (W/mK)'};
if SaveModels == true && optimize_MCMC == 1
    save(fullfile(outDIR, [fname(1:end-4) '_Models.mat']),'models','names')
end
fid = fopen(fullfile(outDIR, fname), 'wt');
if fid == -1
  error('Cannot open log file.');
end
fprintf(fid,'1D thermal model results considering one location (note: %s).\nDatetime: %s Runtime: %0.1f s\n',MData.notes,fname(13:end-4),elapsedTime);
fprintf(fid,'Time step (s): %0.0f\nNsteps: %0.0f\nNwalkers: %0.0f\n# of Layers: %0.0f, Top Layer Size (m):  %0.4f\nInitial Inputs:\n',MData.dt,MData.nstep,MData.nwalkers,length(MData.layer_size),MData.layer_size(1));
for i = 1:MData.nvars
    fprintf(fid,'%s = %0.4f\n',names{:,i},MData.vars_init(i));
end
fprintf(fid,'Goodness of fit: %0.2f\n',fval);
fprintf(fid,'RESULTS [16th pctl, 50th pctl, 84th pctl]:\n');
for i = 1:MData.nvars
    fprintf(fid,'%s = %0.4f, %0.4f, %0.4f\n',names{:,i},RESULTS(:,i));
end
fprintf(fid,'TI = %0.4f, %0.4f, %0.4f\n',TIm,TI,TIp);
fclose(fid);
end