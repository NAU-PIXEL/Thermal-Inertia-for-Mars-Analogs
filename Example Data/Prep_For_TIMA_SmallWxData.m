%Prep_For_TIMA_SmallWxData.m
% Ari Koeppel -- Copyright 2021
%Script to prepare inputs for heat transfer modeling.


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
[LogFile, LogFilePath] = uigetfile('X:\akoeppel\TI_Modeling\TI_Earth_2D\IcelandEmergingEsker2022_10cm\EmergingEsker2022_LoggerData.mat','Select data file containing micrometeorological and thermal data'); %## prompt for log
% ***************

% ***************
%Load data
Data_orig = struct2cell(load(fullfile(LogFilePath, LogFile)));
Data_orig = Data_orig{:,:};
variableNames = Data_orig.Properties.VariableNames;
%Prompt user to select a variable from the table
[selectedIndex, ok] = listdlg('ListString', variableNames, ...
                               'SelectionMode', 'single', ...
                               'PromptString', 'Select a variable:', ...
                               'Name', 'Select Variable');
%Check if the user made a selection
if ok
    selectedVariableName = variableNames{selectedIndex};  % Get the name of the selected variable
    if strcmp('LWLowerCo_Avg',selectedVariableName) || strcmp('LWLowerCo',selectedVariableName)
        Orig_Temps_to_Fit = (Data_orig.LWLowerCo_Avg./5.670374419e-8).^0.25-273.15; %(Data.LWLowerCo_Avg./5.670374419e-8).^0.25-273.15
        Orig_err = 5.*ones(size(Data_orig.LWLowerCo_Avg));
    else
        Orig_Temps_to_Fit = Data_orig.(selectedVariableName);    % Get the values of the selected variable
        Orig_err = (0.002.*(273.15+Orig_Temps_to_Fit)+0.20); %Accuracy is 0.05K after calibration; CS240 Accuracy ± (0.15 + 0.002T)K TO ADD: DIF in temp between back of CS240 and top
    end    
    % Display the selected variable name and its values
    fprintf('You selected: %s\n', selectedVariableName);
else
    fprintf('No variable was selected.\n');
end

observed_TT = timetable(Data_orig.TIMESTAMP,Orig_Temps_to_Fit);
mask = isnan(Orig_Temps_to_Fit);


%% Measured by Hand
[AddlInputs, AIPath] = uigetfile('X:\akoeppel\TI_Modeling\TI_Earth_2D\IcelandEmergingEsker2022_10cm\Additional_Inputs.txt','Select txt file containing density, emissivity, material and max depth inputs');
fileID = fopen(fullfile(AIPath,AddlInputs), 'r');
if fileID == -1
    error('Could not open file: %s', filename);
end
var=1;
while ~feof(fileID)
    line = fgetl(fileID); % Read a line from the file
    if ischar(line)
        % Split the line by the colon
        parts = strsplit(line, ':');
        
        % Check if we got exactly two parts
        if length(parts) == 2
            varValue = strtrim(parts{2});   % Get the value (trim whitespace)
            
            % Convert the value to the appropriate data type
            numericValue = str2double(varValue);
            if ~isnan(numericValue)
                Addldata{var} = numericValue; % Store as a numeric value
                var = var+1;
            else
                Addldata{var} = varValue; % If it's not a number, store as a string
                var = var+1;
            end
        end
    end
end
fclose(fileID);
density = Addldata{1}; %kg/m^3
emissivity = Addldata{2};
T_std = Addldata{3};
material = Addldata{4};%(options: basalt,amorphous,granite,sandstone,clay,salt,ice)
material_lower = Addldata{5};%(options: basalt,amorphous,granite,sandstone,clay,salt,ice)
if ~ismember(material, {'basalt', 'amorphous', 'granite', 'sandstone', 'clay', 'salt', 'ice'}) || ~ismember(material_lower, {'basalt', 'amorphous', 'granite', 'sandstone', 'clay', 'salt', 'ice'})
    error('Error: The material type "%s" is not valid. Please choose from: %s', ...
          materialType, strjoin(validMaterials, ', '));
end
Depth_Max = Addldata{6};%m
evap_depth = Addldata{7}; %Static evap depth (0.075 is reasonable for unsaturated sandy soils)
Vars_init = [Addldata{8};Addldata{9};Addldata{10};Addldata{11};Addldata{12};Addldata{13};Addldata{14};Addldata{15}];
mantle_thickness = Addldata{16}; %minutes
k_dry_std_mantle = Addldata{17};
t_step = Addldata{18}; %minutes
NDAYS = Addldata{19};
nwalkers = Addldata{20}; %100
nstep = Addldata{21}; %10000
burnin_mcmc = Addldata{22}; %fraction of results to omit
burnin_fit = Addldata{23}; %fraction of results to omit
sigma = Addldata{24}; % dictates sinsitivity of walkers to change
nvars = Addldata{25};
minit = zeros(length(Vars_init),nwalkers);
rng(49)  % For reproducibility
for i = 1:nwalkers
    minit(:,i) = Vars_init + sigma*Vars_init.*randn(length(Vars_init),1);
end
VWC_dug_depth = [Addldata{26},Addldata{27},Addldata{28},Addldata{29}]./100; %Depth in m of probe elements

%% Resample
%Error on the FLIR measurements is +/- 5 C or 5% of readings in the -25°C to +135°C range
if t_step >= 1
    Data = retime(Data_orig,'regular','mean','TimeStep',minutes(t_step),'EndValues',NaN);
    if strcmp('LWLowerCo_Avg',selectedVariableName) || strcmp('LWLowerCo',selectedVariableName)
        Temps_to_Fit = (Data.LWLowerCo_Avg./5.670374419e-8).^0.25-273.15; %(Data.LWLowerCo_Avg./5.670374419e-8).^0.25-273.15
        erf = @(T) 5.*ones(length(T));
        err = erf(Data.LWLowerCo_Avg);
    else
        Temps_to_fit = Data.(selectedVariableName);    % Get the values of the selected variable
        erf = @(T) (0.002.*(273.15+T)+0.20);
        err = erf(Temps_to_fit); %Accuracy is 0.05K after calibration; CS240 Accuracy ± (0.15 + 0.002T)K TO ADD: DIF in temp between back of CS240 and top% err = (0.05.*(Temps_to_fit)); %5% or 5K before;
    end
    fit_ind = find(~isnan(Temps_to_fit));

else
    Data = retime(Data_orig,'regular','pchip','TimeStep',minutes(t_step));
    if strcmp('LWLowerCo_Avg',selectedVariableName) || strcmp('LWLowerCo',selectedVariableName)
        Temps_to_Fit = (Data.LWLowerCo_Avg./5.670374419e-8).^0.25-273.15; %(Data.LWLowerCo_Avg./5.670374419e-8).^0.25-273.15
        erf = @(T) 5.*ones(length(T));
        err = erf(Data.LWLowerCo_Avg);
    else
        Temps_to_fit = Data.(selectedVariableName);    % Get the values of the selected variable
        erf = @(T) (0.002.*(273.15+T)+0.20);
        err = erf(Temps_to_fit); %Accuracy is 0.05K after calibration; CS240 Accuracy ± (0.15 + 0.002T)K TO ADD: DIF in temp between back of CS240 and top% err = (0.05.*(Temps_to_fit)); %5% or 5K before;
    end
    observed_TT = timetable(Data_orig.TIMESTAMP,Orig_Temps_to_Fit);
    mask = isnan(Orig_Temps_to_Fit);
    observed_TT(mask,:) = [];
    for times = 1:size(observed_TT,1)
        [~,fit_ind(times)] = min(abs(datenum(Data.TIMESTAMP)-datenum(observed_TT.Time(times))));
    end
end
%Error on the FLIR measurements is +/- 5 C or 5% of readings in the -25°C to +135°C range
testvarnames = {'k-upper' 'Pore network con. par. (mk)' 'Surf. ex. coef. (CH)' 'Surf. ex. coef. (CE)' 'Soil Moist. Infl. (thetak) (%)' 'Soil Moist. Infl. (thetaE) (%)', 'Transition Depth (m)'};
StartTemp_1 = 273.15+Temps_to_fit(1); %Use first observed temperature as start for model top layer
% ***************
%Timestep
dt = seconds(Data.TIMESTAMP(2)-Data.TIMESTAMP(1));
Smooth_Window = 25/(dt/60); %Calculate a reasonable smoothing window for jumpy data (i.e. wind, VWC)
% ***************
%% ***************
TT = timetable(Data.TIMESTAMP,Temps_to_fit);
Interpolated_Temp = fillmissing(TT,'linear');
Interpolated_Temp=Interpolated_Temp{:,1};

% ***************
% Observational Data:
%To Do: Use only SWUpper and individ pixel albedo rather than homogenous SWNet for whole scene.
% R_Short_Net = Data.SWUpper_Avg-Data.SWLower_Avg;    
R_Short_Lower = Data.SWLower_Avg;
R_Short_Lower(R_Short_Lower<0)=0;
R_Short_Upper = Data.SWUpper_Avg;
R_Short_Upper(R_Short_Upper<0)=0;
Albedo = median(Data.SWLower_Avg(Data.SolarElevationCorrectedForAtmRefractiondeg>0)./Data.SWUpper_Avg(Data.SolarElevationCorrectedForAtmRefractiondeg>0)).*ones(size(Data.Albedo_Avg));
R_Long_Upper = Data.LWUpperCo_Avg; %Temp Corrected sky IR radiance
if ~ismember('AirTC', Data.Properties.VariableNames);Data.AirTC = Data.AirTC_Avg;end
Air_Temp_C = Data.AirTC; %Air Temp
Humidity = Data.RH./100; %Rel humid
Pressure_air_Pa = Data.BP_mbar.*100; %barometric pressure - converted to Pa later on
Soil_Temp_C_Probe = Data.T_Avg; %In Situ near surf Soil Temperature from nearest Dry probe
Dewpoint_C = (243.04.*log(Humidity)+17.625.*Air_Temp_C./(243.04+Air_Temp_C))./(17.625-log(Humidity./100)+17.625.*Air_Temp_C./(243.04+Air_Temp_C));
SolarAzimuthCwfromS = 180 - Data.SolarAzimuthAngledegCwFromN;
f_diff = Data.DF;
SolarZenith_Apparent = 90-Data.SolarElevationCorrectedForAtmRefractiondeg;
clear Dug_VWC Dug_Temp
Dug_VWC(:,1) = Data.VWC_Avg;
Dug_VWC(:,2) = Data.VWC_2_Avg;
Dug_VWC(:,3) = Data.VWC_3_Avg;
Dug_Temp(:,1) = Data.T_Avg; %In Situ near surf Soil Temperature from nearest probe
Dug_Temp(:,2) = Data.T_2_Avg; %In Situ near surf Soil Temperature from nearest probe
Dug_Temp(:,3) = Data.T_3_Avg; %In Situ near surf Soil Temperature from nearest probe
if ismember('VWC_4_Avg',Data.Properties.VariableNames)
    Dug_VWC(:,4) = Data.VWC_4_Avg;
    Dug_Temp(:,4) = Data.T_4_Avg; %In Situ near surf Soil Temperature from nearest probe
end
VWC_Smooth_Window = 50/(dt/60);
Dug_VWC_smooth = smoothdata(Dug_VWC,'gaussian',VWC_Smooth_Window);

if ~ismember('WS_ms_Avg', Data.Properties.VariableNames);Data.WS_ms_Avg = Data.mean_wind_spd_ms;end
WindSpeed_ms_10 = Data.WS_ms_Avg; % Wind Speed from sensor @ ~10 ft (3 m) height
WindSpeed_ms_10_smooth = smoothdata(WindSpeed_ms_10,'gaussian',Smooth_Window); %wind is noisy, so this smoothes it

%% Routine to set up boundary and initial conditions for model
FLAY = 0.0075;%0:0.005:1+0.0025*t_step;%:0.0005:0.1;%0:0.005:1+0.0025*t_step;%logspace(-2,0,400);%0:0.005:1; %FLAY values to test Christian used 0.01

% ***************
clear T_Test
%Define Lower Boundary conditions
T_Deep = mean(Dug_Temp(:,end))+273.15;%mean([max(Soil_Temp_C_Dry),min(Soil_Temp_C_Dry)])+273.15; %Static mean of near surf temp

density_plus = density + 997*mean(Dug_VWC(:,end),'omitnan'); %Dry density + water content
Cp_Deep = tima_specific_heat_model_hillel(density,density_plus);
% ***************
RLAY = 1.3;%linspace(1.5,1,0.5/0.05+1); %1.15 Thickness geometric multiplier of layers beneath Christian used 1.3
% ***************
stability_array = FLAY;
check = 0.*[stability_array 1];
check(1)=1; %option to make consecutively stable or not
%Loop to set up grid of layer thicknesses and match each layer to nearest VWC reading
% DSD = sqrt(86400/(pi)*Vars_init(1)/(density*Cp_Deep)); %diurnal Skin depth (meters)
for i = 1:length(stability_array) %Loop to optimize layer thickness (i.e. highest resolution w/o becomming unstable or costing too much computation time)
    % *************** Subsurface Resolution Parameters and Lower Boundary ***************
    %FLAY = fspace(i); %Thickness of emitting top-most layer
    % ***************
    %Loop to set up grid of layer thicknesses and match each layer to nearest VWC reading
    Layer_size_B = FLAY(i); % Top layer thickness
    LayerMult = RLAY;
    if Layer_size_B(1) > 0.10
        error('Layer 1 reached too big at over 10 cm')
    end

    while sum(Layer_size_B) < Depth_Max
        Layer_size_B = cat(2,Layer_size_B,max(Layer_size_B)*LayerMult); %(meters)
    end
    VWC_column = NaN([length(WindSpeed_ms_10) length(Layer_size_B)]);
    for t = 1:length(WindSpeed_ms_10)
        VWC_column(t,:) = interp1(VWC_dug_depth,Dug_VWC_smooth(t,:),cumsum(Layer_size_B),'linear','extrap');
        VWC_column(t,cumsum(Layer_size_B)<=min(VWC_dug_depth)) = Dug_VWC_smooth(t,1);
        VWC_column(t,cumsum(Layer_size_B)>=max(VWC_dug_depth)) = Dug_VWC_smooth(t,end);
    end

    % Initialize Temperatures
    clear Subsurface_Temperatures_Running TEMP T_Start Subsurface_Temperatures
    Subsurface_Temperatures = tima_initialize(Vars_init(1),density,Vars_init(2),Vars_init(5),T_std,T_Deep,Interpolated_Temp,dt,Layer_size_B,VWC_column,Humidity,NDAYS,material,'depth_transition',Vars_init(7),'material_lower',material_lower);
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
        VWC_column,evap_depth.*ones(size(Air_Temp_C)),Humidity,emissivity,...
        Pressure_air_Pa,'albedo',Albedo,'material',material,...
        'depth_transition',FitVar(7),'material_lower',material_lower,'mantle_thickness',mantle_thickness,'k_dry_std_mantle',k_dry_std_mantle);

    % %Need to make this fxn adjustable to dif modes (have it decide
    % %automatically based on which inputs exist?)
    % formod = @(FitVar) tima_heat_transfer(FitVar(1),FitVar(2),FitVar(3),...
    %     FitVar(4),FitVar(5),FitVar(6),density,dt,T_std,Air_Temp_C,R_Short_Upper,...
    %     R_Short_Lower,R_Long_Upper,WindSpeed_ms_10,T_Deep,T_Start,Layer_size_B,...
    %     VWC_column,evap_depth.*ones(size(Air_Temp_C)),Humidity,emissivity,...
    %     Pressure_air_Pa,'albedo',Albedo,'material',material,'depth_transition',FitVar(7),'k_dry_std_lower',FitVar(8),'material_lower',material_lower);
    % %If the initialization looks good, but the MCMC run returns: 'Starting
    % %points for all walkers must have finite logP' then check that the
    % %formod in both scripts is the same.

    Test_Result = formod(Vars_init(:));
    if isreal(Test_Result) && abs(max(Test_Result)-max(Temps_to_fit)) < 50 && sum(sum(max(Subsurface_Temperatures) > (273.15+max(Temps_to_fit(1:(1440/(dt/60))))))) == 0 && sum(isnan(Test_Result)) == 0 %indicators of stability
        check(i+1) = 1;
    end
    if check(i+1) == 1 && (check(i) == check (i+1)) %make sure its consecutively stable
        break
    end
    if i == length(stability_array) && i ~= 1
        error('Max Layer grid size reached with no convergence!')
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
SS= plot(Soil_Temp_C_Probe(1:1440/(dt/60)),'k','LineWidth', 1 ,'DisplayName','near Surf Probe Temp Observed'); %1st day only (1440s)
hold off
legend([M(1:2) SS],'Interpreter','none')
xlabel('Minutes after start time')
ylabel('Temperature (K)')
ttl = sprintf('Day 1 Repeated, %0.f layers, RLAY = %0.2f, FLAY = %0.3f',length(Layer_size_B),LayerMult,FLAY);
title(ttl)
T_Test(:) = T_Start-273.15;%SUB_TEMPS((end-Time_from_end)/(dt/60),:);
figure
plot(cumsum(Layer_size_B)-Layer_size_B./2, T_Test);
hold on
plot(VWC_dug_depth, Dug_Temp((1440)/(dt/60),:));
xlabel('Depth (m)')
ylabel('Temperature (C)')

%% TEST model to make sure it's in the right ball park
% Test_Result = formod(Vars_init(:));
% if sum(isnan(Test_Result))>1 || ~isreal(Test_Result)
%     error('Complex result likely because division by zero in latent heat model')
% end
figure
hold on
M(1) = fill([Data.TIMESTAMP; flipud(Data.TIMESTAMP)],[Temps_to_fit-err;flipud(Temps_to_fit+err)], [128 193 219]./255,'Linestyle','none','DisplayName','FLIR error');
set(M(1), 'edgecolor', 'none');
set(M(1), 'FaceAlpha', 0.5);
G = plot(Data.TIMESTAMP,Soil_Temp_C_Probe,'b','DisplayName','Measured');
H = plot(Data.TIMESTAMP,(Data.LWLowerCo_Avg./5.670374419e-8).^0.25-273.15,'g','DisplayName','CNR4 LWLower, emmis=1');
F = scatter(Data_orig.TIMESTAMP,Orig_Temps_to_Fit,1,'k.','DisplayName','FLIR Observations');
leg = sprintf('Top %0.1f cm modeled temperature', FLAY*100);
M(2)= plot(Data.TIMESTAMP,Test_Result(:,1),'r', 'LineWidth', 2 ,'DisplayName',leg);
hold off
xlabel('Time (hr)');
ylabel('Surface Temperature (C)');
legend([G,H,F,M(1),M(2)],'Interpreter','none')
chi_v_test = sum((Temps_to_fit(fit_ind)-Test_Result(fit_ind)).^2./err(fit_ind).^2)/(length(Temps_to_fit(fit_ind))-nvars); %reduced chi squared
ttl = sprintf('TI top = %0.2f Jm^{-2}K^{-1}s^{-1/2}, chi^{2} = %0.2f', sqrt(Vars_init(1)*Cp_Deep*density),chi_v_test);%Calculate TI from results
title(ttl,'Interpreter','tex','FontName','Ariel')

%% Inputs:
%   Time data - struct of timeseries data variables 
%   Time data - struct of timeseries data variables 
      TData.air_Temp_C=Air_Temp_C;
      TData.DF=f_diff;
      TData.VWC_column=VWC_column;
      TData.evap_depth=evap_depth.*ones(size(Air_Temp_C));
      TData.humidity=Humidity;
      TData.pressure_air_Pa=Pressure_air_Pa;
      TData.r_long_upper=R_Long_Upper;
      TData.r_short_upper=R_Short_Upper;
      TData.r_short_lower=R_Short_Lower;
      TData.solarazimuth_cwfromS=SolarAzimuthCwfromS;
      TData.solarzenith_apparent=SolarZenith_Apparent;
      TData.timed_albedo=Albedo;
      TData.TIMESTAMP = Data.TIMESTAMP;
      TData.temps_to_fit_interp=Interpolated_Temp;
      TData.windspeed_horiz_ms=WindSpeed_ms_10;
      TData.temp_column=Dug_Temp;

%   Model Data - Struct of static and model format variables
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

% clearvars -except out_DIR TData MData
c = fix(clock);
[fname, outDIR] = uiputfile(sprintf('_%02.0f%02.0f-%s.mat',c(4),c(5),date), 'Please navigate to outDIR and enter filename with .mat suffix : ');
save([outDIR,fname],'TData','MData','outDIR','-mat')