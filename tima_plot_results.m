function [] = tima_plot_results(TData,MData,models,names,varargin)
p = inputParser;
p.addRequired('TData',@isstruct);
p.addRequired('MData',@isstruct);
p.addRequired('models');
p.addRequired('names',@ischar);
p.addParameter('TwoSpot',false,@islogical);
p.parse(TData, MData, models, names,varargin{:});
p=p.Results;
TwoSpot = p.TwoSpot;

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
formod_fluxes = @(FitVar) tima_heat_transfer_energy_terms(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.dug_VWC_smooth,TData.evap_depth, MData.VWC_depth_indices,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,'albedo',TData.timed_albedo);
[Result_Temp,Result_latent,Result_keff,q_conv,q_rad,q_G] = formod_fluxes(RESULTS(2,:));
for time = 1:length(TData.temps_to_fit)
    full_latent(time) = sum(Result_latent(time,:));
end
if TwoSpot == true
    formod_fluxes_II = @(FitVar) tima_heat_transfer_energy_terms(FitVar(1),FitVar(2),FitVar(3),...
            FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
            TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
            TData.dug_VWC_smooth_II,TData.evap_depth_II, MData.VWC_depth_indices,TData.humidity,MData.emissivity,...
            TData.pressure_air_Pa,'T_adj1',[2495,300.26],'T_adj2',[2521,301.95],'albedo',TData.timed_albedo_II);
    [Result_Temp_II,Result_latent_II,Result_keff_II,q_conv_II,q_rad_II,q_G_II] = formod_fluxes_II(RESULTS(2,:));
    for time = 1:length(TData.temps_to_fit)
        full_latent_II(time) = sum(Result_latent_II(time,:));
    end
end

figure
hold on
xlabel('Time (hr)');
ylabel('Temperature (C)');
F(1) = fill([TData.TIMESTAMP(MData.fit_ind); flipud(TData.TIMESTAMP(MData.fit_ind))],[TData.temps_to_fit(MData.fit_ind)-TData.err(MData.fit_ind);flipud(TData.temps_to_fit(MData.fit_ind)+TData.err(MData.fit_ind))],[128 193 219]./255,'Linestyle','none','DisplayName','FLIR error');
set(F(1), 'edgecolor', 'none');
set(F(1), 'FaceAlpha', 0.5);
F(2) = scatter(TData.TIMESTAMP(MData.fit_ind),TData.temps_to_fit(MData.fit_ind),1,'k.','DisplayName','FLIR Surface Observations');
M = plot(TData.TIMESTAMP(MData.fit_ind),Result_Temp(MData.fit_ind,1),'r', 'LineWidth', 2 ,'DisplayName','Surface Modeled');

hold off
legend([F(12) M], 'Interpreter','none')
if TwoSpot == true
    chi_v = sum(([TData.temps_to_fit(MData.fit_ind)-tima_formod_subset(RESULTS(2,:),MData.fit_ind,formod); TData.temps_to_fit_II(MData.fit_ind)-tima_formod_subset(RESULTS(2,:),MData.fit_ind,formod_II)]).^2./([TData.err(MData.fit_ind); TData.err_II(MData.fit_ind)]).^2)/(2*length(TData.temps_to_fit(MData.fit_ind))-length(MData.nvars));
else
    chi_v = sum((TData.temps_to_fit(MData.fit_ind)-tima_formod_subset(RESULTS(2,:),MData.fit_ind,formod)).^2./TData.err(MData.fit_ind).^2)/(length(TData.temps_to_fit(MData.fit_ind))-length(MData.nvars));
end
Cp_std = tima_specific_heat_model_hillel(MData.density,MData.density,0);
TI =  sqrt(RESULTS(2,1)*MData.density*Cp_std);
TIp = sqrt(RESULTS(3,1)*MData.density*Cp_std);
TIm = sqrt(RESULTS(1,1)*MData.density*Cp_std);
ttl = sprintf('TI Top [Jm^{-2}K^{-1}s^{-12}] = %0.2f, chi_v = %0.2f',TI,chi_v);%Calculate TI from results
title(ttl,'Interpreter','tex','FontName','Ariel')

% residuals = Temps_to_fit-Test_Result;
% figure(5)
% plot(Data.TIMESTAMP,residuals)
% title('residuals')

%Chi^2_v reduced chi squared = goodness of fit
%-chi^2_v  1 Fit is poor. The chosen model is unlikely to fit the data.
%-chi^2_v  1 Fit does not fully capture the data. It's possible there are too few data points or the uncertainties are underestimated.
%-chi^2_v  1 Fit is too good (overfitting). The model fits noise in the data or the uncertainties are overestimated (can be caused by having too many free parameters in fitting).
%In general chi^2_v ~= 1 means that the fit is doing a good job without overfitting. While these guidelines don't replace full statistical tests, they can be quick check of goodness of fit.

%% Plot sub fluxes!
figure
hold on
ylabel('Temperature (C)');
plot(TData.TIMESTAMP(MData.fit_ind),TData.temps_to_fit(MData.fit_ind),'k','LineWidth', 0.5,'DisplayName','FLIR Surface Observations');
title('FLIR Surface Observations')


figure
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP,q_conv,'g', 'LineWidth', 1,'DisplayName','sensible heat');

figure
hold on
ylabel('Temperature (C)');
plot(TData.TIMESTAMP,Result_keff(:,1),'r', 'LineWidth', 1 ,'DisplayName','k_eff');
hold off

figure
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP,q_rad,'m', 'LineWidth', 1 ,'DisplayName','Radiative heat');

figure
for time = 1:length(TData.temps_to_fit)
    full_latent(time) = sum(Result_latent(time,:));
end
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP,full_latent,'c', 'LineWidth', 1 ,'DisplayName','latent heat');

figure
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP,q_G,'b', 'LineWidth', 1,'DisplayName','ground heat');

end