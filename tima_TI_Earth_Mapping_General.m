function [] = tima_TI_Earth_Mapping_General(TData,MData,inDIR,outDIR,row,varargin)
%% TIMA_TI_EARTH_MAPPING_GENERAL (TI Mars Analogs)
%   Surface energy balance model for deriving thermal inertia in terrestrial sediments using diurnal
%   observations taken in the field to fit 1D multi-parameter model to each pixel. Justification for
%   approach: https://www.mathworks.com/help/gads/table-for-choosing-a-solver.html
%
% Description
%   This model uses a surrogate optimization solver with 1-7 paramters to fit IR-image-observed
%   surface temperatures using meteorological data and surface parameters derived from
%   aerial/satellite imagery. Fitting Parameters can include any or all of: top layer thermal
%   conductivity at 300K, bottome layer thermal conductivity at 300K, Depth of thermal conductivity
%   transition, pore-network-connectivity, Surf. ex. coef. (sensible), Surf. ex. coef. (latent), Soil Moist. Inflection (conductivity) (%), Soil Moist. Infl. (latent heat) (%)
%
% Input Parameters
%   TData: Time data - struct of timeseries data variables 
%       TData.air_Temp_C: 
%       TData.DF:
%       TData.dug_VWC_smooth: 
%       TData.dug_VWC_smooth_II: 
%       TData.evap_depth 
%       TData.evap_depth_II 
%       TData.humidity: as a fraction
%       TData.pressure_air_Pa: 
%       TData.r_long_upper:
%       TData.r_short_upper: 
%       TData.r_short_lower: 
%       TData.solarazimuth_cwfromS:
%       TData.solarzenith_apparent:
%       TData.timed_albedo: 
%       TData.timed_albedo_II: 
%       TData.TIMESTAMP: 
%       TData.temps_to_fit: 
%       TData.temps_to_fit_II: 
%       TData.windspeed_horiz_ms: 
%
%   
%   MData: Model Data - Struct of static and model format variables
%       MData.burnin: 
%       MData.dt:
%       MData.density: 
%       MData.emissivity: 
%       MData.layer_size: 
%       MData.mccount:  NwalkersxNsteps This is the total number of MC proposals, -NOT the number per chain.
%       MData.minit:
%       MData.nvars: 
%       Mdata.parallel: t or f   [f]
%       MData.T_deep: 
%       MData.T_start: 
%       MData.T_std: 
%       MData.ThinChain: [20]
%       MData.VWC_depth_indices: 
%       MData.notes:
%       MData.vars_init:

%   out_DIR:

% Outputs:
%   TK_Line_%u.txt
%   Depth_Line_%u.txt
%   fval_Line_%u.txt
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

format shortG
p = inputParser;
p.addRequired('TData',@isstruct);
p.addRequired('MData',@isstruct);
p.addRequired('outDIR',@ischar);
p.addParameter('Mode','1layer',@ischar);%'2layer','2layer_fixed_depth','2layer_fixed_lower'
p.parse(TData, MData, outDIR, varargin{:});
p=p.Results;
Mode = p.Mode;
%% SET SolarZenith to Apparent (90-apparent solar elevation angle)
%Inputs loaded in .sh file
%load EE2022_Workspace_1Min.mat
format shortG
imptopts = detectImportOptions([inDIR,'Slope.csv']);
imptopts.DataLines = [row row];
Data_Slope_X = readmatrix([inDIR,'Slope.csv'],imptopts);
Data_Aspect_X = readmatrix([inDIR,'Aspect_cwfromS.csv'],imptopts);
Data_Albedo_X = readmatrix([inDIR,'Albedo.csv'],imptopts); 
ShadowDataDir = [inDIR,'Shadows/'];
ShadowFiles = dir(fullfile(ShadowDataDir,'*.csv')); %gets all files
Data_Shadows_X = NaN([size(Data_Albedo_X,2) length(ShadowFiles)]); 
for k = 1:length(ShadowFiles)
    FileName = fullfile(ShadowDataDir, ShadowFiles(k).name);
    Data_Shadows_X(:,k) = readmatrix(FileName,imptopts); 
end
Data_UAV_X = NaN([size(Data_Albedo_X,2) size(MData.UAV_flight_ind,2)]);
for t = 1:size(MData.UAV_flight_ind,2)
    Data_UAV_X(:,t) = readmatrix([inDIR,sprintf('TempC_%u.csv',t)],imptopts);
end
% opts.UseParallel = false;
poolobj = gcp('nocreate');
delete(poolobj);
parpool('Processes',10)
%% Run Model in 1 point mode for each pixel
%row = defined in .sh file
%row = 500
parfor col = 1:MData.col_max
        RESULTS = [NaN NaN];
        fval = NaN;
        if isnan(Data_UAV_X(col,1)) || isnan(Data_UAV_X(col,2)) || isnan(Data_UAV_X(col,3)) ||...
            isnan(Data_UAV_X(col,4)) || isnan(Data_UAV_X(col,5)) || isnan(Data_UAV_X(col,6)) ||...
            isnan(Data_UAV_X(col,7)) || isnan(Data_UAV_X(col,8)) || isnan(Data_UAV_X(col,9)) ||...            
            isnan(Data_Slope_X(col)) || isnan(Data_Aspect_X(col)) || isnan(Data_Albedo_X(col)) || isnan(single(Data_Shadows_X(col,1)))
            continue
        end
        if Mode == '1layer'
            formod = @(theta) tima_heat_transfer(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,MData.vars_init(7),MData.vars_init(8),'albedo',Data_Albedo_X(col),...
                'slope_angle',Data_Slope_X(col),'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',MData.shadow_time_ind,'MappingMode',true);
        elseif Mode == '2layer'
            formod = @(theta) tima_heat_transfer(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,theta(2),theta(3),'albedo',Data_Albedo_X(col),...
                'slope_angle',Data_Slope_X(col),'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',MData.shadow_time_ind,'MappingMode',true,'rho_dry_lower',MData.density_lower);
        elseif Mode == '2layer_fixed_depth'
            formod = @(theta) tima_heat_transfer(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,theta(2),MData.vars_init(8),'albedo',Data_Albedo_X(col),...
                'slope_angle',Data_Slope_X(col),'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',MData.shadow_time_ind,'MappingMode',true,'rho_dry_lower',MData.density_lower);
        elseif Mode == '2layer_fixed_lower'
            formod = @(theta) tima_heat_transfer(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,MData.vars_init(7),theta(2),'albedo',Data_Albedo_X(col),...
                'slope_angle',Data_Slope_X(col),'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',MData.shadow_time_ind,'MappingMode',true,'rho_dry_lower',MData.density_lower);
        end
        Temps_Obs = Data_UAV_X(col,:);
        Temps_Obs = Temps_Obs(:);
        err = (ones(length(Temps_Obs),1)+0.05); 
        Obj = @(theta) sum((Temps_Obs-tima_formod_subset(theta,MData.UAV_flight_ind,formod)).^2./err.^2)/(length(Temps_Obs)-MData.nvars); %Reduced Chi_v         
        problem = struct('solver','surrogateopt','lb',MData.lbound,'ub',MData.ubound,'objective',Obj,'options',MData.probopts,'PlotFcn',[]) ; 
        [RESULTS_holder,fval_holder] = surrogateopt(problem);
        if isempty(RESULTS_holder) || isempty(fval_holder)
            continue
        else
            RESULTS(:) = RESULTS_holder';
            fval = fval_holder;
            writematrix(round(RESULTS(1),3),[outDIR,sprintf('TK_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            writematrix(round(fval,3),[outDIR,sprintf('fval_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            if Mode == '2layer'
                writematrix(round(RESULTS(2),3),[outDIR,sprintf('TK_lower_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
                writematrix(round(RESULTS(3),3),[outDIR,sprintf('Depth_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            elseif Mode == '2layer_fixed_depth'
                writematrix(round(RESULTS(2),3),[outDIR,sprintf('TK_lower_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            elseif Mode == '2layer_fixed_lower'
                writematrix(round(RESULTS(2),3),[outDIR,sprintf('Depth_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            end
        end
end
poolobj = gcp('nocreate');
delete(poolobj);
%% Combine data into rows
if Mode == '1layer'
    LineTKData = NaN([1 MData.col_max]);
    LinefvalData = NaN([1 MData.col_max]);
elseif Mode == '2layer'
    LineTKData = NaN([1 MData.col_max]);
    LineTKlowerData = NaN([1 MData.col_max]);
    LineDepthData = NaN([1 MData.col_max]);
    LinefvalData = NaN([1 MData.col_max]);
elseif Mode == '2layer_fixed_depth'
    LineTKData = NaN([1 MData.col_max]);
    LineTKlowerData = NaN([1 MData.col_max]);
    LinefvalData = NaN([1 MData.col_max]);
elseif Mode == '2layer_fixed_lower'
    LineTKData = NaN([1 MData.col_max]);
    LineDepthData = NaN([1 MData.col_max]);
    LinefvalData = NaN([1 MData.col_max]);
end
for col = 1:MData.col_max
    LineTKFile = [out_DIR,sprintf('TK_Row_%u_Col_%u.txt',row,col)];
    LineTKlowerFile = [out_DIR,sprintf('TK_lower_Row_%u_Col_%u.txt',row,col)];
    LineDepthFile = [out_DIR,sprintf('Depth_Row_%u_Col_%u.txt',row,col)];
    LinefvalFile = [out_DIR,sprintf('fval_Row_%u_Col_%u.txt',row,col)];
    if isfile(LineTKFile)
        LineTKData(1,col) = readmatrix(LineTKFile);
        LinefvalData(1,col) = readmatrix(LinefvalFile);
    end
    if isfile(LineTKlowerFile)
        LineTKlowerData(1,col) = readmatrix(LineTKlowerFile);
    end
    if isfile(LineDepthData)
	    LineDepthData(1,col) = readmatrix(LineDepthFile);
    end
end
writematrix(LineTKData,[outDIR,sprintf('TK_Line_%u.txt',row)],'Delimiter',',')
writematrix(LinefvalData,[outDIR,sprintf('fval_Line_%u.txt',row)],'Delimiter',',')
if Mode == '2layer'
    writematrix(LineTKlowerData,[outDIR,sprintf('TK_lower_Line_%u.txt',row)],'Delimiter',',')
    writematrix(LineDepthData,[outDIR,sprintf('Depth_Line_%u.txt',row)],'Delimiter',',')
elseif Mode == '2layer_fixed_depth'
    writematrix(LineTKlowerData,[outDIR,sprintf('TK_lower_Line_%u.txt',row)],'Delimiter',',')
elseif Mode == '2layer_fixed_lower'
    writematrix(LineDepthData,[outDIR,sprintf('Depth_Line_%u.txt',row)],'Delimiter',',')
end
delete([outDIR,sprintf('*Row_%u_Col*.txt',row)])
end