function [] = tima_TI_Earth_Mapping_General(TData,HData,MData,in_DIR,out_DIR,row)
%% TIMA_TI_EARTH_MAPPING_GENERAL (TI Mars Analogs)
%   Surface energy balance model for deriving thermal inertia in terrestrial sediments using diurnal
%   observations taken in the field to fit 2D multi-parameter model to each pixel. Justification for
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
%       TData.GHI
%       TData.humidity: as a fraction.
%       TData.pressure_air_Pa:
%       TData.r_long_upper:
%       TData.solarazimuth_cwfromS:
%       TData.solarzenith_apparent:
%       TData.windspeed_horiz_ms:
%
%   HData: Hand data - Struct of single measurment variables
%       HData.density:
%       HData.emissivity:
%       HData.reflectivity: average surface reflectivity of site
%       HData.top_start_Temp_K
%   
%   MData: Model Data - Struct of model format variables
%       MData.col_max:
%       MData.dt:
%       MData.layer_size:
%       MData.lbound:
%       MData.nvars:
%       MData.probopts:
%       MData.shadow_time_ind
%       MData.start_Temps_K:
%       MData.UAV_flight_ind:
%       MData.ubound:
%       MData.vars_assigned:
%       MData.VWC_depth_indices:
%       
%
%   in_DIR: directory where raster data is stored (char) - must include Slope.csv, Aspect_cwfromS.csv, Albedo.csv, TempC_1cm_%u.csv, Shadow_Row_%u.csv
%   out_DIR:
%   row: row of raster to be processed


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



%% SET SolarZenith to Apparent (90-apparent solar elevation angle)
%Inputs loaded in .sh file
%load EE2022_Workspace_1Min.mat
format shortG
imptopts = detectImportOptions([DIR,'Slope_1cm.csv']);
imptopts.DataLines = [row row];
Data_Slope_X = readmatrix([DIR,'Slope_1cm.csv'],imptopts);
Data_Aspect_X = readmatrix([DIR,'Aspect_cwfromS_1cm.csv'],imptopts);
Data_Albedo_X = readmatrix([DIR,'Albedo_1cm.csv'],imptopts); %1x8840
Data_Shadow_X = readmatrix([DIR,'Shadows_1cm_rows/',sprintf('Shadow_Row_%u.csv',row)]);
Data_UAV_X = NaN([size(Data_Albedo_X,2) size(MData.UAV_flight_ind,2)]);
for t = 1:size(MData.UAV_flight_ind,2)
    Data_UAV_X(:,t) = readmatrix([DIR,sprintf('TempC_1cm_%u.csv',t)],imptopts);
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
            isnan(Data_Slope_X(col)) || isnan(Data_Aspect_X(col)) || isnan(Data_Albedo_X(col)) || isnan(single(Data_Shadow_X(col,1)))
            continue
        end
        formod = @(theta) Heat_Transfer(theta(1),MData.vars_assigned(2),theta(2),MData.vars_assigned(4),MData.vars_assigned(5),MData.vars_assigned(6),MData.vars_assigned(7),MData.vars_assigned(8),HData.density,MData.dt,TData.air_Temp_C,TData.GHI,TData.DF,TData.r_long_upper,...
            TData.windspeed_horiz_ms,HData.top_start_Temp_K,MData.start_Temps_K,MData.layer_size,TData.dug_VWC_smooth,MData.VWC_depth_indices,TData.humidity,HData.emissivity,Data_Albedo_X(col),Data_Slope_X(col),Data_Aspect_X(col),...
            TData.solarzenith_apparent,TData.solarazimuth_cwfromS,TData.pressure_air_Pa,single(Data_Shadow_X(col,:)),MData.shadow_time_ind,HData.reflectivity);%%
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
            writematrix(round(RESULTS(1),3),[OUT_DIR,sprintf('TK_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            writematrix(round(RESULTS(2),3),[OUT_DIR,sprintf('Depth_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            writematrix(round(fval,3),[OUT_DIR,sprintf('fval_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
        end
end

%% Combine data into rows

LineTKData = NaN([1 MData.col_max]);
LineDepthData = NaN([1 MData.col_max]);
LinefvalData = NaN([1 MData.col_max]);
for col = 1:MData.col_max
    LineTKFile = [OUT_DIR,sprintf('TK_Row_%u_Col_%u.txt',row,col)];
    LineDepthFile = [OUT_DIR,sprintf('Depth_Row_%u_Col_%u.txt',row,col)];
    LinefvalFile = [OUT_DIR,sprintf('fval_Row_%u_Col_%u.txt',row,col)];
    if isfile(LineTKFile)
        LineTKData(1,col) = readmatrix(LineTKFile);
	    LineDepthData(1,col) = readmatrix(LineDepthFile);
        LinefvalData(1,col) = readmatrix(LinefvalFile);
    end
end
writematrix(LineTKData,[OUT_DIR,sprintf('TK_Line_%u.txt',row)],'Delimiter',',')
writematrix(LineDepthData,[OUT_DIR,sprintf('Depth_Line_%u.txt',row)],'Delimiter',',')
writematrix(LinefvalData,[OUT_DIR,sprintf('fval_Line_%u.txt',row)],'Delimiter',',')
delete([OUT_DIR,sprintf('*Row_%u_Col*.txt',row)])
poolobj = gcp('nocreate');
delete(poolobj);
end