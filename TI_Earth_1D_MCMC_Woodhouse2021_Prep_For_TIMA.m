%TI_Earth_1D_MCMC.m
% Ari Koeppel -- Copyright 2021
%Script to model thermal inertia in terrestrial sediments using diurnal 
%observations taken in the field fit to multi-parameter model using MCMC 
%simulations. This model applies to one point/pixel (1D).

%Justification for approach: https://www.mathworks.com/help/gads/table-for-choosing-a-solver.html
%SURROGATE OPTIMIZATION IS USEFUL FOR TESTING
%MCMC IS USEFUL FOR OBTAINING UNCERTAINTY ON TI/k results

%References: Kieffer, 2013, 

%% Clear Workspace
clc;
clear
close all;
format compact; format longG;

%% Load Measured Inputs
% ***************
%Logger data
%Find and load a .mat file with data already synchronized and
%clipped to the proper time window. See script: PSTAR_Data_Integration_SunsetMay2021_Moving.m
% Retrieve data file path name
LogFile = dir('X:\common\FIELD_CAMPAIGNS\WoodhouseMesa_May2021\Ground_Station\Synced_Data\WoodhouseMesa_May2021_Irradiance.mat');
% ***************

% ***************
%Load data
Data_orig = struct2cell(load(fullfile(LogFile.folder, LogFile.name)));
Data_orig = Data_orig{:,:};
Orig_Temps_to_Fit = Data_orig.FLIR_DrySurf_Corr;
Orig_err = (0.002.*(273.15+Orig_Temps_to_Fit)+0.20); %Accuracy is 0.05K after calibration; CS240 Accuracy ± (0.15 + 0.002T)K TO ADD: DIF in temp between back of CS240 and top
observed_TT = timetable(Data_orig.TIMESTAMP,Orig_Temps_to_Fit);
mask = isnan(Orig_Temps_to_Fit);
starts = [mask(1); (diff(mask)>0)];
stops = [(diff(mask)<0);~mask(end)];
t_step = 0.25; %minutes
if t_step >= 1
    Data = retime(Data_orig,'regular','mean','TimeStep',minutes(t_step));
else
    Data = retime(Data_orig,'regular','pchip','TimeStep',minutes(t_step));
end
Data.mask = zeros(height(Data),1);
Data.mask(Data_orig.TIMESTAMP(starts & ~stops)) = 1; % leave singleton blocks alone
Data.mask(Data_orig.TIMESTAMP(stops & ~starts)) = -1;% leave singleton blocks alone
Data.mask = cumsum(Data.mask);
Data.mask(Data_orig.TIMESTAMP(stops)) = 1; % mark singleton blocks as missing
% ind = find(Data.mask==0);
observed_TT(mask,:) = [];
for times = 1:size(observed_TT,1)
    [~,fit_ind(times)] = min(abs(datenum(Data.TIMESTAMP)-datenum(observed_TT.Time(times))));
end
% ***************

% ***************
%Timestep
dt = seconds(Data.TIMESTAMP(2)-Data.TIMESTAMP(1));
Smooth_Window = 25/(dt/60); %Calculate a reasonable smoothing window for jumpy data (i.e. wind, VWC)
% ***************

% ***************
%Measured by Hand
k_meas_avg = 0.115; % SA1: 0.106, SA2: 0.115
% TI_meas_avg = SA1: 149.1 SA2: 307.5; %Jm^{-2}K^{-1}s^{-1/2}
density = 1487; %kg/m^3 SA1.1: 1232.16, SA1.2: 1487.23, SA2: 1351.33
emissivity = 0.966; %SA_Avg:0.966, SA_Basalt: 0.961, SA_Sandstone: 0.944
Med_Grain = 0.002; %m from sieving
porosity = 1-density/2900;
k_H2O = 0.61;
T_std = 300; %K
material = 'basalt';
MappingMode = 0;
NDAYS = 20;
Depth_Max = 0.5;
% ***************

%% Extracting Variables
% ***************
% ROIs:
Corrected_Dry_Point_Temp = Data.FLIR_DrySurf_Corr;
Corrected_Wet_Point_Temp = Data.FLIR_WetSurf_Corr;
% ***************
TT = timetable(Data.TIMESTAMP,Corrected_Dry_Point_Temp);
Interpolated_Temp = fillmissing(TT,'linear');
Interpolated_Temp=Interpolated_Temp{:,1};

% ***************
% Observational Data:
%To Do: Use only SWUpper and individ pixel albedo rather than homogenous SWNet for whole scene.
R_Short_Net = Data.SWUpper_Avg_AllData-Data.SWLower_Avg_AllData;    
R_Short_Lower = Data.SWLower_Avg_AllData;%-min(Data.SWLower_Avg_AllData);
R_Short_Lower(R_Short_Lower<0)=0;
R_Short_Upper = Data.SWUpper_Avg_AllData;%-min(Data.SWUpper_Avg_AllData);
R_Short_Upper(R_Short_Upper<0)=0;
Albedo = median(Data.SWLower_Avg_AllData(Data.SolarElevationCorrectedForAtmRefractiondeg>0)./Data.SWUpper_Avg_AllData(Data.SolarElevationCorrectedForAtmRefractiondeg>0));%Data.Albedo_Avg;
% Albedo = median(Data.SWLower_Avg_AllData./Data.SWUpper_Avg_AllData);%Data.Albedo_Avg_AllData;
R_Long_Upper = Data.LWUpperCo_Avg_AllData; %Temp Corrected sky IR radiance
Air_Temp_C = Data.AirTC; %Air Temp
Humidity = Data.RH./100; %Rel humid
Pressure_air_Pa = Data.BP_mbar.*100; %barometric pressure - converted to Pa later on
Soil_Temp_C_Dry = Data.SoilTemp2; %In Situ near surf Soil Temperature from nearest Dry probe
Soil_Temp_C_Wet = Data.SoilTemp6; %In Situ near surf Soil Temperature from nearest Wet probe
Dewpoint_C = (243.04.*log(Humidity./100)+17.625.*Air_Temp_C./(243.04+Air_Temp_C))./(17.625-log(Humidity./100)+17.625.*Air_Temp_C./(243.04+Air_Temp_C));
SolarAzimuthCwfromS = 180 - Data.SolarAzimuthAngledegCwFromN;
f_diff = Data.DF;
SolarZenith_Apparent = 90-Data.SolarElevationCorrectedForAtmRefractiondeg;
Aspect_cwfromS = 0;


VWC_Smooth_Window = 100/(dt/60);
Dug_VWC_Wet(:,1) = Data.VWC6; %In Situ near surf Soil Moisture from nearest WET probe
Dug_VWC_Wet(:,2) = Data.VWC_101;
Dug_VWC_Wet(:,3) = Data.VWC_201;
Dug_VWC_Wet(:,4) = Data.VWC_301;
Dug_VWC_Wet(:,5) = Data.VWC_401;
Dug_VWC_Wet(:,6) = Data.VWC_501;
Dug_VWC_Wet_smooth = smoothdata(Dug_VWC_Wet,'gaussian',VWC_Smooth_Window);
SoilVWC_Wet_smooth = smoothdata(Data.VWC6,'gaussian',Smooth_Window); %VWC is discretized, so this smoothes it
Dug_VWC_Wet_smooth(:,1) = SoilVWC_Wet_smooth;

Dug_VWC_Dry(:,1) = Data.VWC_56; %In Situ near surf Soil Moisture from nearest DRY probe
Dug_VWC_Dry(:,2) = Data.VWC_106;
Dug_VWC_Dry(:,3) = Data.VWC_203;
Dug_VWC_Dry(:,4) = Data.VWC_303;
Dug_VWC_Dry(:,5) = Data.VWC_403;
Dug_VWC_Dry(:,6) = Data.VWC_503;
Dug_VWC_Dry_smooth = smoothdata(Dug_VWC_Dry,'gaussian',VWC_Smooth_Window);
% SoilVWC_Dry_smooth = smoothdata(Data.VWC2,'gaussian',Smooth_Window); %VWC is discretized, so this smoothes it
% Dug_VWC_Dry_smooth(:,2) = SoilVWC_Dry_smooth;

% if Surface == 1
TF_char = 'Data.FLIR_WetandDrySurf_Corr'
    Temps_to_fit_II=Data.FLIR_WetSurf_Corr; %Define Wet fitting data here
    Timed_Albedo_II = Data.FLIR_VIS_albedo_wet;
    Soil_Temp_II = Soil_Temp_C_Wet;
    Dug_VWC_smooth_II = Dug_VWC_Wet_smooth;
    Dug_VWC_II = Dug_VWC_Wet;
% else
    Temps_to_fit=Data.FLIR_DrySurf_Corr; %Define Wet fitting data here
    Timed_Albedo = Data.FLIR_VIS_albedo_dry;
    Soil_Temp = Soil_Temp_C_Dry;
    Dug_VWC_smooth = Dug_VWC_Dry_smooth;
    Dug_VWC = Dug_VWC_Dry;
% end

MaxVWC = max(max(Dug_VWC_Wet)); %Saturation at site of wetting;
theta_E = 0.75*MaxVWC;%max(Data.VWC1); %Saturation at site of wetting;
WindSpeed_ms_10 = Data.WS_ms_10; % Wind Speed from sensor @ ~10 ft (3 m) height
% WindSpeed_ms_10_smooth = smoothdata(WindSpeed_ms_10,'gaussian',Smooth_Window); %wind is noisy, so this smoothes it

WindSpeed_ms_30 = Data.WS_ms_30; % Wind Speed from sensor @ ~20 ft (6 m) height
WindSpeed_ms_30_smooth = smoothdata(WindSpeed_ms_30,'gaussian',Smooth_Window); %wind is noisy, so this smoothes it

% WindSpeed_deviation = WindSpeed_ms_30_smooth-WindSpeed_ms_10_smooth-mean(WindSpeed_ms_30_smooth-WindSpeed_ms_10_smooth,'omitnan');

% Buried VWC probe data from nearest probe
% VWC_dug_depth = [5,10,20,30,40,50]; %Depth in cm of probe elements
VWC_dug_depth = [0.01,5,15,25,35,45]; %Depth in cm of probe elements


%Error on the FLIR measurements is +/- 5 C or 5% of readings in the -25°C to +135°C range
err = (0.002.*(273.15+Temps_to_fit)+0.20); %Accuracy is 0.05K after calibration; CS240 Accuracy ± (0.15 + 0.002T)K TO ADD: DIF in temp between back of CS240 and top
% err = (0.05.*(Temps_to_fit)); %5% or 5K before;
err_II = (0.002.*(273.15+Temps_to_fit_II)+0.20); %Accuracy is 0.05K after calibration; CS240 Accuracy ± (0.15 + 0.002T)K
N_VWC = length(VWC_dug_depth);
evap_depth_wet = ones(size(Temps_to_fit));
evap_depth_dry = ones(size(Temps_to_fit));
for t = 1:length(Temps_to_fit)
    for z = 1:N_VWC-1
        dVWC(z) = (Dug_VWC_Wet_smooth(t,z+1)-Dug_VWC_Wet_smooth(t,z))./(VWC_dug_depth(z+1)-VWC_dug_depth(z));
        if z == 1 && Dug_VWC_Wet_smooth(t,1) > 0 && dVWC(1)<= 0
            evap_depth_wet(t) = VWC_dug_depth(1);
            break
        elseif z > 1 && dVWC(z)-dVWC(z-1) < 0
            evap_depth_wet(t) = VWC_dug_depth(z);
            break
        end
    end
    for z = 1:N_VWC-1
        dVWC(z) = (Dug_VWC_Dry_smooth(t,z+1)-Dug_VWC_Dry_smooth(t,z))./(VWC_dug_depth(z+1)-VWC_dug_depth(z));
        if z == 1 && Dug_VWC_Dry_smooth(t,1) > 0 && dVWC(1)<= 0
            evap_depth_dry(t) = VWC_dug_depth(1);
            break
        elseif z > 1 && dVWC(z)-dVWC(z-1) < 0
            evap_depth_dry(t) = VWC_dug_depth(z);
            break
        end
    end
end
%% Test Variables [k above transition depth; k below transition depth; Sensible Heat Multiplier];
Vars_init = [0.17;0.73;410;3100;0.35;0.06];%[0.15;0.23;0.992;420;2700;0.66;0.07];
names = {'k-upper''Pore network con. par. (mk)' 'Surf. ex. coef. (CH)' 'Surf. ex. coef. (CE)' 'Soil Moist. Infl. (thetak) (%)' 'Soil Moist. Infl. (thetaE) (%)'};
StartTemp_1 = 273.15+Temps_to_fit(1); %Use first observed temperature as start for model top layer

%% Routine to set up boundary and initial conditions for model
fspace = 0.005+0.0025*t_step;%logspace(-2,0,400);%0:0.005:1; %FLAY values to test Christian used 0.01
check = 0.*fspace;
% ***************
%Define Lower Boundary conditions
T_Deep = mean([max(Soil_Temp_C_Dry),min(Soil_Temp_C_Dry)])+273.15; %Static mean of near surf temp

density_plus = density + 997*mean(Dug_VWC_smooth(:,end),'omitnan'); %Dry density + water content
Cp_Deep = tima_specific_heat_model_hillel(density,density_plus,mean(Dug_VWC_smooth(:,end),'omitnan'));
% ***************
RLAY = 1.3; %1.15 Thickness geometric multiplier of layers beneath Christian used 1.3
% ***************
%Loop to set up grid of layer thicknesses and match each layer to nearest VWC reading
DSD = sqrt(86400/(pi)*Vars_init(1)/(density*Cp_Deep)); %diurnal Skin depth (meters)
for i = 1:length(fspace) %Loop to optimize layer thickness (i.e. highest resolution w/o becomming unstable or costing too much computation time)
    % *************** Subsurface Resolution Parameters and Lower Boundary ***************
    FLAY = fspace(i); %Thickness of emitting top-most layer
    % ***************
    %Loop to set up grid of layer thicknesses and match each layer to nearest VWC reading
    Layer_size_B = FLAY; % Top layer thickness
    % if Layer_size_B(1) < Med_Grain*2
    %     continue
    % end
    if Layer_size_B(1) > 0.10
        error('Layer 1 reached too big at over 10 cm')
    end
    count = 2;
    [M,VWC_depth_indices(1)] = min(abs(Layer_size_B-VWC_dug_depth./100));
    while sum(Layer_size_B) < Depth_Max
      Layer_size_B = cat(2,Layer_size_B,max(Layer_size_B)*RLAY); %(meters)
      [M,VWC_depth_indices(count)] = min(abs(sum(Layer_size_B)-VWC_dug_depth./100));
      count = count+1;
    end
    % ***************
    % Initialize Temperatures
    clear Subsurface_Temperatures_Running TEMP T_Start Subsurface_Temperatures
    Subsurface_Temperatures = tima_initialize(Vars_init(1),density,Vars_init(2),Vars_init(5),T_std,T_Deep,Interpolated_Temp,dt,Layer_size_B,Dug_VWC_Dry_smooth,VWC_depth_indices,Humidity,NDAYS,material);
    Subsurface_Temperatures_Running(:,:) = Subsurface_Temperatures(1,:,:)-273.15; %Set up array using first day
    for D = 2:size(Subsurface_Temperatures,1) % for plotting
        TEMP(:,:) = Subsurface_Temperatures(D,:,:)-273.15;
        Subsurface_Temperatures_Running = cat(1,Subsurface_Temperatures_Running,TEMP);
    end
    T_Start(:)=Subsurface_Temperatures(end,end,:);
    T_Start(1) = StartTemp_1;

    formod = @(FitVar) tima_heat_transfer(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),density,dt,T_std,Air_Temp_C,R_Short_Upper,...
        R_Short_Lower,R_Long_Upper,WindSpeed_ms_10,T_Deep,T_Start,Layer_size_B,...
        Dug_VWC_smooth,evap_depth_dry,VWC_depth_indices,Humidity,emissivity,...
        Pressure_air_Pa,'albedo',Timed_Albedo);

    formod_II = @(FitVar) tima_heat_transfer(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),density,dt,T_std,Air_Temp_C,R_Short_Upper,...
        R_Short_Lower,R_Long_Upper,WindSpeed_ms_10,T_Deep,T_Start,Layer_size_B,...
        Dug_VWC_smooth_II,evap_depth_wet,VWC_depth_indices,Humidity,emissivity,...
        Pressure_air_Pa,'T_adj1',[2495,300.26],'T_adj2',[2521,301.95],'albedo',Timed_Albedo_II);

    Test_Result = formod(Vars_init(:));
    if isreal(Test_Result) && abs(max(Test_Result)-max(Temps_to_fit)) < 50 && sum(sum(max(Subsurface_Temperatures) > (273.15+max(Temps_to_fit(1:(1440/(dt/60))))))) == 0 && sum(isnan(Test_Result)) == 0 %indicators of stability
        check(i+1) = 1;
    end
    if check(i) == 1 && (check(i) == check (i+1)) %make sure its consecutively stable
        break
    end
end

%% Plot Initialization Results to see if converges
figure
M(1)= plot(Subsurface_Temperatures_Running((end-1440/(dt/60)+1):end,1),'r', 'LineWidth', 2 ,'DisplayName','Surface Observed');
hold on
for i = 2:size(Subsurface_Temperatures_Running,2)
    M(i) = plot(Subsurface_Temperatures_Running((end-1440/(dt/60)+1):end,i),'b--', 'LineWidth', 2,'DisplayName','Subsurface_Modeled');
hold on
end
SS= plot(Soil_Temp_C_Dry(1:1440/(dt/60)),'k','LineWidth', 1 ,'DisplayName','Depth Temp Observed'); %1st day only (1440s)
hold off
legend([M(1:2) SS],'Interpreter','none')
xlabel('Minutes after start time')
ylabel('Temperature (K)')
ttl = sprintf('Day 1 Repeated, %0.f layers, RLAY = %0.2f, FLAY = %0.2f',length(Layer_size_B),RLAY,FLAY);
title(ttl)
T_Test(:) = T_Start-273.15;%SUB_TEMPS((end-Time_from_end)/(dt/60),:);
figure
plot([0 cumsum(Layer_size_B)], [T_Test T_Test(end)]);
xlabel('Depth (m)')
ylabel('Temperature (C)')

% %% Input observations into heat transfer function to generate forward model
% formod = @(theta) Heat_Transfer(theta(1),theta(2),theta(3),theta(4),theta(5),theta(6),density,dt,Air_Temp_C,R_Short_Net,R_Long_Upper,...
%     WindSpeed_ms_30_smooth,T_Deep,StartTemp_1,T_Start,Layer_size_B,Dug_VWC_smooth,VWC_depth_indices,Humidity./100,emissivity,Timed_Albedo,Pressure_air.*100);%

%% TEST model to make sure it's in the right ball park
% Test_Result = formod(Vars_init(:));
if sum(isnan(Test_Result))>1 || ~isreal(Test_Result)
    error('Complex result likely because division by zero in latent heat model')
end
figure
hold on
M(1) = fill([Data.TIMESTAMP(fit_ind); flipud(Data.TIMESTAMP(fit_ind))],[Orig_Temps_to_Fit(~mask)-Orig_err(~mask);flipud(Orig_Temps_to_Fit(~mask)+Orig_err(~mask))], [128 193 219]./255,'Linestyle','none','DisplayName','FLIR error');
set(M(1), 'edgecolor', 'none');
set(M(1), 'FaceAlpha', 0.5);
G = plot(Data.TIMESTAMP,Soil_Temp,'b','DisplayName','Measured');
F = scatter(Data_orig.TIMESTAMP,Orig_Temps_to_Fit,1,'k.','DisplayName','FLIR Observations');
leg = sprintf('Top %0.1f cm modeled temperature', FLAY*DSD*100);
M(2)= plot(Data.TIMESTAMP,Test_Result(:,1),'r', 'LineWidth', 2 ,'DisplayName',leg);
hold off
xlabel('Time (hr)');
ylabel('Surface Temperature (C)');
legend([G,F,M(1),M(2)],'Interpreter','none')
chi_v_test = sum((Temps_to_fit(fit_ind)-Test_Result(fit_ind)).^2./err(fit_ind).^2)/(length(Temps_to_fit(fit_ind))-length(Vars_init)); %reduced chi squared
ttl = sprintf('TI top = %0.2f Jm^{-2}K^{-1}s^{-1/2}, chi^{2} = %0.2f', sqrt(Vars_init(1)*Cp_Deep*density),chi_v_test);%Calculate TI from results
title(ttl,'Interpreter','tex','FontName','Ariel')

%% Initialize Walkers
% ***************
%For the old and the new model parameters, determine
% the posterior (the log likelihood + the log prior) and find
% R = Pnew/Pold (which is lnR = lnPnew - lnPold in log space)
% – If R > 1 always accept the new parameters
% – If R < 1 accept the new parameters with probability R
% – Check how often the new parameters are accepted. If
% this is far from ~30% (meaning inefficient), change the proposal step size (sigma)
% Simulation Parameters
nwalkers = 10; %100
nstep = 100; %10000
mccount = nstep*nwalkers;% This is the total number, -NOT the number per chain.% What is the desired total number of monte carlo proposals.
burnin = 0.5; %fraction of results to omit
sigma = 10^-3; % dictates sinsitivity of walkers to change
rng(49)  % For reproducibility
minit = zeros(length(Vars_init),nwalkers);
for i = 1:nwalkers
    minit(:,i) = Vars_init + sigma*Vars_init.*randn(length(Vars_init),1);
end
nvars = 6;
%% Inputs:
%   Time data - struct of timeseries data variables 
      TData.air_Temp_C=Air_Temp_C;
      TData.DF=f_diff;
      TData.dug_VWC_smooth=Dug_VWC_smooth;
      TData.dug_VWC_smooth_II=Dug_VWC_smooth_II;
      TData.err = err;
      TData.err_II = err_II;
      TData.evap_depth=evap_depth_dry;
      TData.evap_depth_II =evap_depth_wet;
      TData.humidity=Humidity;
      TData.pressure_air_Pa=Pressure_air_Pa;
      TData.r_long_upper=R_Long_Upper;
      TData.r_short_upper=R_Short_Upper;
      TData.r_short_lower=R_Short_Lower;
      TData.solarazimuth_cwfromS=SolarAzimuthCwfromS;
      TData.solarzenith_apparent=SolarZenith_Apparent;
      TData.timed_albedo=Timed_Albedo;
      TData.timed_albedo_II=Timed_Albedo_II;
      TData.TIMESTAMP = Data.TIMESTAMP;
      TData.temps_to_fit=Temps_to_fit;
      TData.temps_to_fit_II=Temps_to_fit_II;
      TData.windspeed_horiz_ms=WindSpeed_ms_10;

%   Model Data - Struct of static and model format variables
      MData.burnin=burnin;
      MData.dt=dt;  
      MData.density=density; 
      MData.emissivity=emissivity;
      MData.fit_ind = fit_ind;
      MData.layer_size=Layer_size_B;
      MData.mccount = mccount;
      MData.minit = minit;
      MData.nvars=nvars;
      MData.parallel= false;
      MData.T_deep= T_Deep; 
      MData.T_start= T_Start;
      MData.T_std=T_std;
      MData.ThinChain=20; %[20]
      MData.VWC_depth_indices = VWC_depth_indices;
      MData.notes = 'WH2021 Wet & Dry';
      MData.vars_init = Vars_init;

      outDIR='X:\akoeppel\TI_EARTH_1D\Woodhouse2021';

% clearvars -except out_DIR TData MData
c = fix(clock);                       
fname = sprintf('WH2021Inputs_%02.0f%02.0f-%s.mat',c(4),c(5),date);
save([outDIR,'\',fname],'TData','MData','outDIR','-mat')