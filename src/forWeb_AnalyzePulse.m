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
% get values of coefficients
% First-order rate coefficients, k (1/yr)
k_mat = forWeb_rate_coeffs(k_factors, Lriver_FHgP, IHgD_pristine, IHgP_pristine);

% translate k values to correct coefficients
k_A_oHgII = k_mat(1);
k_A_tHgII = k_mat(2);
k_A_tHg0 = k_mat(3);
k_A_oHg0 = k_mat(4);
fdep_tf = k_mat(5);
fdep_ts = k_mat(6);
fdep_ta = k_mat(7);
k_Oc_ev = k_mat(8);
k_Oc_sp1 = k_mat(9);
k_Oc_vsi = k_mat(10);
k_Oc_sp2 = k_mat(11);
k_Oc_vis = k_mat(12);
k_Oc_vid = k_mat(13);
k_Oc_sp3 = k_mat(14);
k_Oc_vdi = k_mat(15);
k_Te_rf = k_mat(16);
k_Te_p = k_mat(17);
k_T_exfs = k_mat(18);
k_T_exfa = k_mat(19);
k_Te_BBf_1 = k_mat(20);
k_Te_BBf_2 = k_mat(21);
k_Te_rs = k_mat(22);
k_T_exsf = k_mat(23);
k_T_exsa = k_mat(24);
k_Te_BBs_1 = k_mat(25);
k_Te_BBs_2 = k_mat(26);
k_Te_ra = k_mat(27);
k_T_exaf = k_mat(28);
k_T_exam = k_mat(29);
k_Te_BBa_1 = k_mat(30);
k_Te_BBa_2 = k_mat(31);
k_T_riv_f = k_mat(32);
k_T_riv_s = k_mat(33);
k_T_riv_a = k_mat(34);
k_O_riv_f = k_mat(35);
k_O_riv_s = k_mat(36);
k_O_riv_a = k_mat(37);
E_geo = k_mat(38);
f_HgPexport = k_mat(39);
% Set biomass burning for anthropogenic era
k_Te_BBf = k_Te_BBf_2; % fast soil
k_Te_BBs = k_Te_BBs_2; % slow soil
k_Te_BBa = k_Te_BBa_2; % armored soil    

% set up time variable
y_end = 2110;
tspan    = -2000:dt:(y_end + 0.8);         % run all the way through end year
t = tspan;
% set pulse size equal to river or atmospheric pulse
if  strcmp(Lpulse,'atmpulse')
    pulse_size_T = pulse_size;
elseif  strcmp(Lpulse,'riverpulse')
    pulse_size_T = river_pulse;
end

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
if Lplot
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
end
%%
% calculate fits for decay
% assume that pulse time is halfway through pulse year, calculate time since pulse for normalization
t_pulse = t - (pulse_time + 0.5);
% take data starting in year following pulse for making fit
indx_fit_start = find(round(t,1) == pulse_time + 1 );
if (strcmp(Lpulse,'riverpulse')) % only have a river pulse, no  atmos pulse
    indx_fit_start = find(round(t,1) == pulse_time + 1.6 ); % ignore first 1.5 year where have jagged response
end
indx_fit_end = find(round(t,1) == pulse_time + 100.8 );

t_fit = t_pulse(indx_fit_start:indx_fit_end);

annual_t = 0.5:100.5;

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

% ouptut coeffs in correct order, since it's arbitrary whether smallest 
% lifetime is first or second
ind_min = find(coeffs_e2 == min(coeffs_e2(2), coeffs_e2(4)));% find smaller exponential lifetime term
ind_max = find(coeffs_e2 == max(coeffs_e2(2), coeffs_e2(4)));% find larger exponential lifetime term

coeffs_emis = [coeffs_e2(ind_min-1), -coeffs_e2(ind_min), coeffs_e2(ind_max-1), -coeffs_e2(ind_max)];
%%
if Lplot
    % make plot
    % colors to use
    c_o = [0.3010 0.7450 0.9330]; % color for ocean
    c_t = [0.4660 0.6740 0.1880]; % color for terrestrial
    c3 = [0.8500, 0.3250, 0.0980];

    % plot
    figure('Position', [100 100 1400 500])
    subplot(1,3,1)
    set(gca, 'FontSize',15)
    hold on
    t_plot_pulse = annual_t + pulse_time + 0.5;% after calculation, shift time so can plot it
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
end