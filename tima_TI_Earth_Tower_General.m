%% TIMA_TI_EARTH_TOWER_GENERAL (TI Mars Analogs)
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
Data_Shadow_X = readmatrix([Shadow_DIR,sprintf('Shadow_Row_%u.csv',row)]);
Data_UAV_X = NaN([size(Data_Albedo_X,2) size(ind,2)]);
for t = 1:size(ind,2)
    Data_UAV_X(:,t) = readmatrix([DIR,sprintf('TempC_1cm_%u.csv',t)],imptopts);
end
% opts.UseParallel = false;
%% Initialize Walkers
% ***************
% Simulation Parameters
nwalkers = 50; %100
nstep = 5000; %10000
burnin = 0.5; %fraction of results to omit
sigma = 10^-2; % dictates sinsitivity of walkers to change
rng(49)  % For reproducibility
minit = zeros(length(Vars_init),nwalkers);
for i = 1:nwalkers
    minit(:,i) = Vars_init + sigma*Vars_init.*randn(length(Vars_init),1);
end
% ***************
%% Inputs to emcee
err_all = [err; err_II];
logPfuns = {@(theta)tima_ln_prior(theta) @(theta)-0.5*sum(([Temps_to_fit-tima_formod_subset(theta,ind,formod); Temps_to_fit_II-tima_formod_subset(theta,ind,formod_II)]).^2./err_all.^2)};% a cell of function handles returning the log probality of a each outcome
mccount = nstep*nwalkers;% What is the desired total number of monte carlo proposals.
%                            This is the total number, -NOT the number per chain.

%% Run emcee
% ***************
% Named Parameter-Value pairs:
%   'StepSize': unit-less stepsize (default=2.5).
%   'ThinChain': Thin all the chains by only storing every N'th step (default=10)
%   'ProgressBar': Show a text progress bar (default=true)
%   'Parallel': Run in ensemble of walkers in parallel. (default=false)
%   'BurnIn': fraction of the chain that should be removed. (default=0)
% ***************
tic
[models,LogPs]=gwmcmc(minit,logPfuns,mccount,'BurnIn',burnin,'Parallel',true,'ThinChain',20);
elapsedTime = toc
if (size(models,1)<size(models,2))&&(ismatrix(models)), models=models'; end %Consider this behaviour further....
if ndims(models)==3
    models=models(:,:)'; 
end
M=size(models,2);
for r=1:M
quant=quantile(models(:,r),[0.16,0.5,0.84]);
RESULTS(:,r) = quant;
end
%%
Result = formod(RESULTS(2,:));
chi_v = sum(([Temps_to_fit-tima_formod_subset(RESULTS(2,:),ind,formod); Temps_to_fit_II-tima_formod_subset(RESULTS(2,:),ind,formod_II)]).^2./err_all.^2)/(2*length(Temps_to_fit)-length(Vars_init)); %reduced chi squared/(length(Temps_to_fit)-length(Vars_init)); %reduced chi squared
TI =  sqrt(RESULTS(2,1)*density*Cp_std);
TIp = sqrt(RESULTS(3,1)*density*Cp_std);
TIm = sqrt(RESULTS(1,1)*density*Cp_std);

%% Print Log File 
    c = fix(clock);
    fname = sprintf('%02.0f%02.0f-%s_VWC_0.txt',c(4),c(5),date);
    fid = fopen(fullfile(LogFile.folder, fname), 'wt');
    if fid == -1
      error('Cannot open log file.');
    end
    fprintf(fid,'Thermal model results of %s from Data: %s\nDatetime: %s Runtime: %0.1f s\n',TF_char,fullfile(LogFile.name),fname(1:end-4),elapsedTime);
    fprintf(fid,'Nsteps: %0.0f\nNwalkers: %0.0f\n# of Layers: %0.0f\nInitial Inputs:\n',nstep,nwalkers,length(Layer_size_B));
    for i = 1:length(Vars_init)
        fprintf(fid,'%s = %0.4f\n',names{:,i},Vars_init(i));
    end
    fprintf(fid,'Reduced Chi Squared: %0.2f\n',chi_v);
    fprintf(fid,'RESULTS [16th pctl, 50th pctl, 84th pctl]:\n');
    for i = 1:length(Vars_init)
        fprintf(fid,'%s = -%0.4f%% %0.4f +%0.4f%%\n',names{:,i},abs(RESULTS(1,i)-RESULTS(2,i))/RESULTS(2,i)*100,RESULTS(2,i),abs(RESULTS(3,i)-RESULTS(2,i))/RESULTS(2,i)*100);
    end
    fprintf(fid,'TI = %0.4f, %0.4f, %0.4f\n',TIm,TI,TIp);
    fclose(fid);
