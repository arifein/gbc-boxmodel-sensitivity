%==========================================================================
% OBJECTIVE
%   Analyze the pulse experiment and calculate coefficients for EAMD and
%   EAME
% REVISION HISTORY
%   20 Jul 2012 - hma - modified from run_1450to2008.m in version 5 of the
%                       code. Take out the mineral reservoir and treat 
%                       geogenic emissions as an external forcing. 
%   26 Jul 2012 - hma - emit anthropogenic emissions at a constant rate of
%                       the course of a single year rather than
%                       interpolating at a sub-annual scale
%   19 Aug 2012 - hma - add dep_contrib as an output
%   06 Feb 2013 - hma - add additional pre-1850 scenarios to do an
%                       uncertainty analysis for anthropogenic emissions
%   06 Feb 2013 - hma - add atmospheric reservoir (Matm) as an output. Add
%                       Lwrite, filepath, and Ldisp as inputs. 
%   30 Jan 2014 - HMA - convert from a function to a standard module, makes
%                       it easier to use. Move AnthroEmiss call to main.m
%   30 Apr 2014 - HMA - make future projections compatible w/ updated box
%                       model version with rivers 
%   08 Sep 2014 - HMA - clean up code and comments for public release
%   28 Mar 2023 - AF - edit code for pulse emissions setup
% Helen M. Amos, hamos@hsph.harvard.edu
%==========================================================================

% set pulse size equal to sum of river and atmospheric pulse (one should be
% zero
pulse_size_T = pulse_size + river_pulse;
disp('River pulse is ')
disp(river_pulse)
disp('Atmos pulse is ')
disp(pulse_size)

% calculate differences in each reservoir after pulse
pulse_diff =  M_pulse(:,:) - M_steady(:,:); 
% Matm = M(1,:); % atmosphere
% Mtf  = M(2,:); % fast soil pool
% Mts  = M(3,:); % slow soil pool
% Mta  = M(4,:); % armored soil pool
% Mocs = M(5,:); % surface ocean
% Moci = M(6,:); % intermediate ocean
% Mocd = M(7,:); % deep ocean

% 1) plot the burden spike and disappearance from atmosphere, ocean, and quick land
figure
indx_pulse = find(round(t,1) == pulse_time-1);
indx_pulse_end = find(round(t,1) == pulse_time+100.8); % arbitrary start date

set(gca, 'FontSize',13)
hold on
plot (t(indx_pulse:indx_pulse_end), pulse_diff(1,indx_pulse:indx_pulse_end), 'k', 'linewidth', 2.5)
title ('Surface Hg Reservoirs - Pulse')
plot (t(indx_pulse:indx_pulse_end), pulse_diff(2,indx_pulse:indx_pulse_end),'linestyle','-.', 'linewidth',2.5,'Color',[0.7 0.7 0.7])
plot (t(indx_pulse:indx_pulse_end), pulse_diff(5,indx_pulse:indx_pulse_end),'linestyle','--','linewidth',2.5,'color',[0.4 0.4 0.4])
legend ('atmosphere', 'fast terrestrial', 'surface ocean', 'Location', 'NorthEast')
xlabel('Time (years)')
ylabel('Mg of Hg')
xlim([t(indx_pulse) t(indx_pulse_end)])
hold off
%%
% 2) plot and calculate EAMD using 2-term exponential, try 2 ways
% calculate deposition
deposition_diff = pulse_diff(1,:) * (k_A_oHgII + k_A_oHg0 + k_A_tHgII + k_A_tHg0);

% calculate fits for decay
% assume that pulse time is halfway through pulse year, calculate time since pulse for normalization
t_pulse = t - (pulse_time + 0.5);
% take data starting in year following pulse for making fit
indx_fit_start = find(round(t,1) == pulse_time + 1 );
if (pulse_size == 0) % only have a river pulse, no  atmos pulse
    indx_fit_start = find(round(t,1) == pulse_time + 1.6 ); % ignore first 1.5 year where have jagged response
end
indx_fit_end = find(round(t,1) == pulse_time + 100.8 );

dep_diff_fit = deposition_diff(indx_fit_start:indx_fit_end) / pulse_size_T; % normalize by pulse size
t_fit = t_pulse(indx_fit_start:indx_fit_end);

% selin version of equation
fo = fitoptions('Method','NonlinearLeastSquares','Lower', [0, -Inf, 0, -Inf], 'Upper', [Inf, 0, Inf, 0]);
ft1 = fittype("a*exp(b*x) + c*(1-exp(d*x))", 'options', fo);
[fit1, gof1] = fit(t_fit(:), dep_diff_fit(:), ft1); 
coeffs1 = coeffvalues(fit1);
% calculate annual values from fit, accounting for pulse size
annual_t = 0.5:100.5;
yfit1 = pulse_size_T * (coeffs1(1) * exp(coeffs1(2) * annual_t) + coeffs1(3) * (1 - exp(coeffs1(4) * annual_t))); 
t_plot_pulse = annual_t + pulse_time + 0.5;% after calculation, shift time so can plot it
strfit1 = "fit v1: R^2 = " + sprintf('%0.3f',gof1.rsquare);
streq1 = sprintf('%0.2f',coeffs1(1)) + "exp(" + sprintf('%0.2f',coeffs1(2)) ...
    + "t) + " + sprintf('%0.2f',coeffs1(3)) + "(1-exp(" + ...
    sprintf('%0.2f',coeffs1(4)) + "t))";

% selin 2018 coefficients
yfit2 = pulse_size_T * (0.68 * exp(-0.43 * annual_t) + 0.008 * (1 - exp(-1.3 * annual_t))); 

% exp2 version of equation
fo = fitoptions('Method','NonlinearLeastSquares','Lower', [0, -Inf, 0, -Inf], 'Upper', [Inf, 0, Inf, 0]);
ft3 = fittype("a*exp(b*x) + c*exp(d*x)", 'options', fo);
[fit3, gof3] = fit(t_fit(:), dep_diff_fit(:), ft3); 
coeffs3 = coeffvalues(fit3);
yfit3 = pulse_size_T * (coeffs3(1) * exp(coeffs3(2) * annual_t) + coeffs3(3) * exp(coeffs3(4) * annual_t)); 
strfit3 = "fit v2: R^2 = " + sprintf('%0.3f',gof3.rsquare);
streq3 = sprintf('%0.2f',coeffs3(1)) + "exp(" + sprintf('%0.2f',coeffs3(2)) ...
    + "t) + " + sprintf('%0.2f',coeffs3(3)) + "exp(" + ...
    sprintf('%0.2f',coeffs3(4)) + "t)";


% make plot
% colors to use
c1 = [0, 0.4470, 0.7410];
c2 = [0.9290, 0.6940, 0.1250];
c3 = [0.8500, 0.3250, 0.0980];
% plot
figure('Position', [100 100 1400 500])
subplot(1,3,1)
set(gca, 'FontSize',15)
hold on
plot (t(indx_pulse:indx_pulse_end), deposition_diff(1,indx_pulse:indx_pulse_end), 'k', 'linewidth', 3)
plot(t_plot_pulse,yfit1,'o','Color',c1, 'MarkerFaceColor',c1)
plot(t_plot_pulse,yfit2,'o','Color',c2, 'MarkerFaceColor',c2)
plot(t_plot_pulse,yfit3,'o','Color',c3, 'MarkerFaceColor',c3)
title ('Atmospheric deposition from pulse')
legend ('Amos model', '2-term exponential fit v1','2-term exponential fit Selin18', '2-term exponential fit v2', 'Location', 'NorthEast')
xlabel('Time (years)')
ylabel('Mg yr^{-1}')
% add text for labels
text(0.15, 0.78, streq1, 'Color',c1, 'Fontsize',15, 'Units','Normalized')
text(0.15, 0.72, strfit1, 'Color',c1, 'Fontsize',15, 'Units','Normalized')
text(0.15, 0.64, streq3, 'Color',c3, 'Fontsize',15, 'Units','Normalized')
text(0.15, 0.58, strfit3, 'Color',c3, 'Fontsize',15, 'Units','Normalized')
xlim([t(indx_pulse) t(indx_pulse_end)])
hold off
% same plot, just zoomed in to the beginning
subplot(1,3,2)
set(gca, 'FontSize',15)
hold on
plot (t(indx_pulse:indx_pulse_end), deposition_diff(1,indx_pulse:indx_pulse_end), 'k', 'linewidth', 3)
plot(t_plot_pulse,yfit1,'o','Color',c1, 'MarkerFaceColor',c1)
plot(t_plot_pulse,yfit2,'o','Color',c2, 'MarkerFaceColor',c2)
plot(t_plot_pulse,yfit3,'o','Color',c3, 'MarkerFaceColor',c3)
title ('Deposition from pulse (zoom begin)')
xlabel('Time (years)')
ylabel('Mg yr^{-1}')
xlim([t(indx_pulse) t(indx_pulse)+21])
hold off

% same plot, just zoomed in to the end
subplot(1,3,3)
set(gca, 'FontSize',15)
hold on
plot (t(indx_pulse:indx_pulse_end), deposition_diff(1,indx_pulse:indx_pulse_end),  'linewidth', 3)
plot(t_plot_pulse,yfit1,'o','Color',c1, 'MarkerFaceColor',c1)
plot(t_plot_pulse,yfit2,'o','Color',c2, 'MarkerFaceColor',c2)
plot(t_plot_pulse,yfit3,'o','Color',c3, 'MarkerFaceColor',c3)
title ('Deposition from pulse (zoom end)')
xlabel('Time (years)')
ylabel('Mg yr^{-1}')
xlim([t(indx_pulse_end)-31 t(indx_pulse_end)])
hold off
%%
% 3) plot and calculate EAME using double exponential, single exponential
% calculate emissions difference in pulse experiment
emiss_diff = pulse_diff(2,:) * (k_Te_rf + k_Te_p + k_Te_BBf) + ... % fast ter
             pulse_diff(3,:) * (k_Te_rs + k_Te_BBs) + ... % slow ter
             pulse_diff(4,:) * (k_Te_ra + k_Te_BBa) + ...% armoured ter
             pulse_diff(5,:) * (k_Oc_ev); % surface ocean
emiss_diff_o = pulse_diff(5,:) * (k_Oc_ev); % ocean contribution
emiss_diff_l = emiss_diff - emiss_diff_o; % land contribution
% take data starting in year following pulse for making fit
emiss_diff_fit = emiss_diff(indx_fit_start:indx_fit_end) / pulse_size_T; % normalize by pulse size

% exp1 version of equation
fo = fitoptions('Method','NonlinearLeastSquares','Lower', [0, -Inf], 'Upper', [Inf, 0]);
ft_e1 = fittype("a*exp(b*x)", 'options', fo);
[fit_e1, gof_e1] = fit(t_fit(:), emiss_diff_fit(:), ft_e1); 
coeffs_e1 = coeffvalues(fit_e1);
yfit_e1 = pulse_size_T * (coeffs_e1(1) * exp(coeffs_e1(2) * annual_t)); 
strfit_e1 = "fit exp1: R^2 = " + sprintf('%0.3f',gof_e1.rsquare);
streq_e1 = sprintf('%0.2f',coeffs_e1(1)) + "exp(" + sprintf('%0.2f',coeffs_e1(2)) ...
    + "t)";

% exp2 version of equation
fo = fitoptions('Method','NonlinearLeastSquares','Lower', [0, -Inf, 0, -Inf], 'Upper', [Inf, 0, Inf, 0]);
ft_e2 = fittype("a*exp(b*x) + c*exp(d*x)", 'options', fo);
[fit_e2, gof_e2] = fit(t_fit(:), emiss_diff_fit(:), ft_e2); 
coeffs_e2 = coeffvalues(fit_e2);
yfit_e2 = pulse_size_T * (coeffs_e2(1) * exp(coeffs_e2(2) * annual_t) + coeffs_e2(3) * exp(coeffs_e2(4) * annual_t)); 
strfit_e2 = "fit v2: R^2 = " + sprintf('%0.3f',gof_e2.rsquare);
streq_e2 = sprintf('%0.2f',coeffs_e2(1)) + "exp(" + sprintf('%0.2f',coeffs_e2(2)) ...
    + "t) + " + sprintf('%0.2f',coeffs_e2(3)) + "exp(" + ...
    sprintf('%0.2f',coeffs_e2(4)) + "t)";

% Selin18 version of equation with 2 exponentials
ft_e2_v1 = fittype("a*exp(b*x) + c*(1-exp(d*x))", 'options', fo);
[fit_e2_v1, gof_e2_v1] = fit(t_fit(:), emiss_diff_fit(:), ft_e2_v1); 
coeffs_e2_v1 = coeffvalues(fit_e2_v1);
% calculate annual values from fit, accounting for pulse size
yfit_e2_v1 = pulse_size_T * (coeffs_e2_v1(1) * exp(coeffs_e2_v1(2) * annual_t) + coeffs_e2_v1(3) * (1 - exp(coeffs_e2_v1(4) * annual_t))); 
strfit_e2_v1 = "fit v1: R^2 = " + sprintf('%0.3f',gof1.rsquare);
streq_e2_v1 = sprintf('%0.2f',coeffs_e2_v1(1)) + "exp(" + sprintf('%0.2f',coeffs_e2_v1(2)) ...
    + "t) + " + sprintf('%0.2f',coeffs_e2_v1(3)) + "(1-exp(" + ...
    sprintf('%0.2f',coeffs_e2_v1(4)) + "t))";

%%
% make plot
% colors to use
c_o = [0.3010 0.7450 0.9330]; % color for ocean
c_t = [0.4660 0.6740 0.1880]; % color for terrestrial
% plot
figure('Position', [100 100 1400 500])
subplot(1,3,1)
set(gca, 'FontSize',15)
hold on
plot (t(indx_pulse:indx_pulse_end), emiss_diff(1,indx_pulse:indx_pulse_end), 'k', 'linewidth', 3)
plot (t(indx_pulse:indx_pulse_end), emiss_diff_o(1,indx_pulse:indx_pulse_end), 'Color', ...
    c_o, 'linewidth', 2)
plot (t(indx_pulse:indx_pulse_end), emiss_diff_l(1,indx_pulse:indx_pulse_end), 'Color', ...
    c_t, 'linewidth', 2)
plot(t_plot_pulse,yfit_e2,'o','Color',c3, 'MarkerFaceColor',c3)

title ('Legacy emissions from pulse')
legend ('Amos model', 'Ocean contribution','Terrestrial contribution', ...
     '2-term exponential fit v2',... 
    'Location', 'NorthEast')
xlabel('Time (years)')
ylabel('Mg yr^{-1}')
% add text for labels
text(0.15, 0.78, streq_e2, 'Color',c3, 'Fontsize',15, 'Units','Normalized')
text(0.15, 0.72, strfit_e2, 'Color',c3, 'Fontsize',15, 'Units','Normalized')
xlim([t(indx_pulse) t(indx_pulse_end)])
hold off
subplot(1,3,2)
set(gca, 'FontSize',15)
hold on
plot (t(indx_pulse:indx_pulse_end), emiss_diff(1,indx_pulse:indx_pulse_end), 'k', 'linewidth', 3)
plot (t(indx_pulse:indx_pulse_end), emiss_diff_o(1,indx_pulse:indx_pulse_end), 'Color', ...
    c_o, 'linewidth', 2)
plot (t(indx_pulse:indx_pulse_end), emiss_diff_l(1,indx_pulse:indx_pulse_end), 'Color', ...
    c_t, 'linewidth', 2)
plot(t_plot_pulse,yfit_e2,'o','Color',c3, 'MarkerFaceColor',c3)

title ('Legacy emissions from pulse (zoom begin)')
xlabel('Time (years)')
ylabel('Mg yr^{-1}')
% add text for labels
xlim([t(indx_pulse) t(indx_pulse)+21])
hold off

subplot(1,3,3)
set(gca, 'FontSize',15)
hold on
plot (t(indx_pulse:indx_pulse_end), emiss_diff(1,indx_pulse:indx_pulse_end), 'k', 'linewidth', 3)
plot (t(indx_pulse:indx_pulse_end), emiss_diff_o(1,indx_pulse:indx_pulse_end), 'Color', ...
    c_o, 'linewidth', 2)
plot (t(indx_pulse:indx_pulse_end), emiss_diff_l(1,indx_pulse:indx_pulse_end), 'Color', ...
    c_t, 'linewidth', 2)
plot(t_plot_pulse,yfit_e2,'o','Color',c3, 'MarkerFaceColor',c3)

title ('Legacy emissions from pulse (zoom end)')
xlabel('Time (years)')
ylabel('Mg yr^{-1}')
% add text for labels
xlim([t(indx_pulse_end)-31 t(indx_pulse_end)])
hold off

disp(fit_e2)