function [] = tima_plot_results(TData,MData,models,names,varargin)
p = inputParser;
p.addRequired('TData',@isstruct);
p.addRequired('MData',@isstruct);
p.addRequired('models');
p.addRequired('names',@iscellstr);
p.addParameter('Wetting',false,@islogical);
p.parse(TData, MData, models, names,varargin{:});
p=p.Results;
Wetting = p.Wetting;

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
formod = @(FitVar) tima_heat_transfer_skin(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,'albedo',TData.timed_albedo,'material',MData.material);

formod_fluxes = @(FitVar) tima_heat_transfer_skin_energy_terms(FitVar(1),FitVar(2),FitVar(3),...
        FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
        TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
        TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
        TData.pressure_air_Pa,'albedo',TData.timed_albedo);
[Result_Temp_Surf,Result_Temp_Sub,q_latent,Result_keff,q_conv,q_rad,q_G] = formod_fluxes(RESULTS(2,:));
if Wetting == true
    formod_fluxes = @(FitVar) tima_heat_transfer_skin_energy_terms(FitVar(1),FitVar(2),FitVar(3),...
            FitVar(4),FitVar(5),FitVar(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
            TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
            TData.VWC_II_column,TData.evap_depth_II,TData.humidity,MData.emissivity,...
            TData.pressure_air_Pa,'T_adj1',[2495,300.26],'T_adj2',[2521,301.95],'albedo',TData.timed_albedo_II);
    [Result_Temp_Surf,Result_Temp_Sub,q_latent,Result_keff,q_conv,q_rad,q_G] = formod_fluxes(RESULTS(2,:));
end

figure
hold on
xlabel('Time (hr)');
ylabel('Temperature (C)');
if Wetting == true
    F(1) = fill([TData.TIMESTAMP(MData.fit_ind); flipud(TData.TIMESTAMP(MData.fit_ind))],[TData.temps_to_fit_II(MData.fit_ind)-TData.err_II(MData.fit_ind);flipud(TData.temps_to_fit_II(MData.fit_ind)+TData.err_II(MData.fit_ind))],[128 193 219]./255,'Linestyle','none','DisplayName','FLIR error');
    F(2) = scatter(TData.TIMESTAMP(MData.fit_ind),TData.temps_to_fit_II(MData.fit_ind),1,'k.','DisplayName','FLIR Surface Observations');
else
    F(1) = fill([TData.TIMESTAMP(MData.fit_ind); flipud(TData.TIMESTAMP(MData.fit_ind))],[TData.temps_to_fit(MData.fit_ind)-TData.err(MData.fit_ind);flipud(TData.temps_to_fit(MData.fit_ind)+TData.err(MData.fit_ind))],[128 193 219]./255,'Linestyle','none','DisplayName','FLIR error');
    F(2) = scatter(TData.TIMESTAMP(MData.fit_ind),TData.temps_to_fit(MData.fit_ind),1,'k.','DisplayName','FLIR Surface Observations');
end
set(F(1), 'edgecolor', 'none');
set(F(1), 'FaceAlpha', 0.5);
M = plot(TData.TIMESTAMP,Result_Temp_Surf(:,1),'r', 'LineWidth', 2 ,'DisplayName','Surface Modeled');

hold off
legend([F(2) M], 'Interpreter','none')
chi_v = sum((TData.temps_to_fit(MData.fit_ind)-tima_formod_subset(RESULTS(2,:),MData.fit_ind,formod)).^2./TData.err(MData.fit_ind).^2)/(length(TData.temps_to_fit(MData.fit_ind))-length(MData.nvars))
Cp_std = tima_specific_heat_model_hillel(MData.density,MData.density);
TI =  sqrt(RESULTS(2,1)*MData.density*Cp_std)
TIp = sqrt(RESULTS(3,1)*MData.density*Cp_std);
TIm = sqrt(RESULTS(1,1)*MData.density*Cp_std);
ttl = sprintf('TI Top [Jm^{-2}K^{-1}s^{-12}] = %0.2f, chi_v = %0.2f',TI,chi_v);%Calculate TI from results
title(ttl)

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
if Wetting == true
    plot(TData.TIMESTAMP(MData.fit_ind),TData.temps_to_fit_II(MData.fit_ind),'k','LineWidth', 0.5,'DisplayName','FLIR Surface Observations');
else
    plot(TData.TIMESTAMP(MData.fit_ind),TData.temps_to_fit(MData.fit_ind),'k','LineWidth', 0.5,'DisplayName','FLIR Surface Observations');
end
title('FLIR Surface Observations')


figure
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP,q_conv,'g', 'LineWidth', 1,'DisplayName','sensible heat');
title('Sensible')

figure
hold on
ylabel('W/mK');
plot(TData.TIMESTAMP,Result_keff(:,1),'r', 'LineWidth', 1 ,'DisplayName','k_eff');
hold off
title('K_eff')

figure
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP,q_rad,'m', 'LineWidth', 1 ,'DisplayName','Radiative heat');
title('Radiative')


figure
for time = 1:length(TData.temps_to_fit)
    full_latent(time) = sum(q_latent(time,:));
end
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP,full_latent,'c', 'LineWidth', 1 ,'DisplayName','latent heat');
% hold on
% plot(TData.TIMESTAMP,q_latent(:,7),'b', 'LineWidth', 0.5 ,'DisplayName','latent heat');
title('Latent')

figure
hold on
ylabel('W/m^2');
plot(TData.TIMESTAMP,q_G,'b', 'LineWidth', 1,'DisplayName','ground heat');
title('Ground')


end