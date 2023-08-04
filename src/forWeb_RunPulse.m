%==========================================================================
% OBJECTIVE
%   Simulate the global perturbation introduced by anthropogenic mercury 
%   sources. Global reservoirs of Hg are initialized from natural
%   stead-state levels, then this module simulates the anthropogenic era. 
%
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

%for safety's sake
clear Matm Mtf Mts Mta Mocs Moci Mocd E1;   


%--------------------------------------------------------------------------
% SET UP
%--------------------------------------------------------------------------

% Assemble matrix A, such that dM/dt = A*M + E
sim_type = 2;
forWeb_makeA

% Emissions (Mg/yr)
E1 = diag(ones(1,7),0);  % indentity matrix
Eg = [E_geo; 0; 0; 0; 0; 0; 0]; % geogenic emissions to the atmosphere

% integration time (yrs)
annual_dt = 1/dt;                   % number of time steps in a year
tspan    = -2000:dt:(y_end + 0.8);         % run all the way through end year

%--------------------------------------------------------------------------
% Hg discharge from rivers
%--------------------------------------------------------------------------
     
% time span and time step you want to interpolate to
t_river_inc = t_SF(1):dt:(2008 + 0.8); % increasing river input
t_river = t_SF(1):dt:(y_end + 0.8);

% intialize
rivHgP_MgYr  = zeros( 1, numel( t_river )); % HgP inputs to each basin (Mg/yr) each decade
rivHgD_MgYr  = zeros( 1, numel( t_river )); % HgD inputs to each basin (Mg/yr) each decade

% interpolate
% global HgD and HgP inputs ocean margins (Mg/yr)
rivHgP_MgYr(1:length(t_river_inc))  = pchip( t_SF, sum(river_HgP_MgYr_save,1), t_river_inc );
rivHgD_MgYr(1:length(t_river_inc))  = pchip( t_SF, sum(river_HgD_MgYr_save,1), t_river_inc );

% set steady after 2008
rivHgP_MgYr(length(t_river_inc)+1:end) = rivHgP_MgYr(length(t_river_inc));
rivHgD_MgYr(length(t_river_inc)+1:end) = rivHgD_MgYr(length(t_river_inc));

% find indices of time where want to input river pulse throughout 2010
pulse_idx_r_s = find(round(t_river,1) == pulse_time ); % start index
pulse_idx_r_e = find(round(t_river,1) == pulse_time + 1 - dt ); % end index
% add river pulse for experiment
%frac_D = 0.043 / 100.; % fraction of released mercury going to dissolved phase
frac_D = 50 / 100.; % fraction of released mercury going to dissolved phase
frac_P = 1 - frac_D; % fraction of released mercury going to particulate phase
rivHgP_MgYr(pulse_idx_r_s:pulse_idx_r_e) = rivHgP_MgYr(pulse_idx_r_s:pulse_idx_r_e) + frac_P * river_pulse; % add release to particulate
rivHgD_MgYr(pulse_idx_r_s:pulse_idx_r_e) = rivHgD_MgYr(pulse_idx_r_s:pulse_idx_r_e) + frac_D * river_pulse; % add release to dissolved

% for storing quasi-direct anthropogenic contribution below when you
% calculate M(t)
store_Mriv_quasi_margin = zeros(1,numel(t_river));

% for storing total riverine discharges to ocean margins below when you
% calculate M(t)
store_Mriv_total_margin      = zeros(1,numel(t_river));
store_Mriv_background_margin = zeros(1,numel(t_river));
store_coastal_burial         = zeros(1,numel(t_river));


%--------------------------------------------------------------------------
% Anthropogenic emissions
%--------------------------------------------------------------------------

% Anthropogenic emissions are interpolated from a decadal to annual scale.
% Emit anthropogenic Hg at a constant rate over a year. 

clear y; % for saftey's sake
AnthroTemp = []; % intialize

for y = 1:length(FTime);
    AnthroTemp            = vertcat(AnthroTemp           , Anthro(y)*ones(annual_dt,1));
end

Anthro            = AnthroTemp;

% hfig = figure(81+ff);
% ff = ff+1;
% set(hfig,'units','normalized','Position',[0.1 0.4 0.5 .7])
% set(gcf,'Color',[1 1 1])
% set(gca,'FontSize',18)
% hold on;
% plot(tspan,Anthro ,'b','LineWidth',3)
% xlabel('Year (AD)')
% ylabel('(Mg a^{-1}) ')
% xlim([-2000 y_end])
% title('Primary Anthropogenic Emissions for Pulse Case')
% hold off;
%%
%--------------------------------------------------------------------------
% Solve M(t) forward in time, stop at end year
%--------------------------------------------------------------------------

% dummy matrix of zeros (not necessary, but dramatically saves time)
M = zeros(7, numel(tspan));

% Initial conditions (Mg)
M(:,1) = [Ratm_PI; Rtf_PI; Rts_PI; Rta_PI; Rocs_PI; Roci_PI; Rocd_PI];

% counter for rivers
jRiv = 1;

for j = 2:numel(tspan); 
    
    % Time depdendent anthropogenic emissions (Mg/yr)
    Ea = [Anthro(j); 0; 0; 0; 0; 0; 0]; 
    
    % M(t + dt) = (A*M(t) + E)*dt + M(t-1)
    M(:,j) = ( A*M(:,j-1) + E1*(Ea+Eg) )*dt + M(:,j-1); 
    
    % for safety's sake
    clear Ea Ep;
    
    % Quasi-direct anthropogenic contribution to rivers (i.e., the
    % difference between total observed riverine discharges and the
    % background contribution).
    
    
    % The quasi-direct anthropogenic contribution comes from products,
    % which didn't become important until 1850. So before 1850 rivers
    % are 100% background
    if tspan(j) >= 1850;
        
        % for safety's sake
        clear Mriv_background_margin  Mriv_background_ocean ...
              Mriv_total_margin       Mriv_total_ocean      ...
              Mriv_quasi_margin       Mriv_quasi_ocean;
        
        % background contribution
        
        % background riverine discharges to ocean margins and open ocean (Mg)
        Mriv_background_margin = dt*(k_T_riv_f*M(2,j) + k_T_riv_s*M(3,j) + k_T_riv_a*M(4,j));
        Mriv_background_ocean  = dt*(k_O_riv_f*M(2,j) + k_O_riv_s*M(3,j) + k_O_riv_a*M(4,j));
        
        % total riverine discharged to ocean margins and reaching the open ocean (Mg)
        Mriv_total_margin      = dt*(rivHgD_MgYr(jRiv) + rivHgP_MgYr(jRiv));
        Mriv_total_ocean       = dt*(rivHgD_MgYr(jRiv) + f_HgPexport*rivHgP_MgYr(jRiv));
        
        % quasi-direct anthropogenic contribution to riverine discharges
        % to ocean margins and what reaches the open ocean (Mg)
        Mriv_quasi_margin      = Mriv_total_margin - Mriv_background_margin;
        Mriv_quasi_ocean       = Mriv_total_ocean  - Mriv_background_ocean;
        
        % prevent quasi-direct anthropogenic contribution from going
        % negative
        if Mriv_quasi_margin < 0;
            Mriv_quasi_margin = 0;
        end
        if Mriv_quasi_ocean < 0;
            Mriv_quasi_ocean  = 0;
        end
        
        % add quasi-direct anthropogenic contribution to open ocean (Mg)
        M(5,j) = M(5,j) + Mriv_quasi_ocean;
        
        % store the quasi-direct contribution to total discharges to
        % ocean margins (Mg), so you can compare to product releases to
        % land and water
        store_Mriv_quasi_margin(jRiv) = Mriv_quasi_margin;
        
        store_Mriv_total_margin(jRiv) = Mriv_total_margin; % total discharges to coastal margins
        store_Mriv_background_margin(jRiv) = Mriv_background_margin; % background discharges to coastal margins
        
        % amount of anthro Hg removed to coastal benthic sediment (Mg)
        store_coastal_burial(jRiv) = Mriv_quasi_margin - Mriv_quasi_ocean;

        % increment counter
        jRiv = 1 + jRiv;
    end
    
end

%%
% parse output
Matm = M(1,:); % atmosphere
Mtf  = M(2,:); % fast soil pool
Mts  = M(3,:); % slow soil pool
Mta  = M(4,:); % armored soil pool
Mocs = M(5,:); % surface ocean
Moci = M(6,:); % intermediate ocean
Mocd = M(7,:); % deep ocean

%--------------------------------------------------------------------------
% PLOTS
%--------------------------------------------------------------------------

% time vector, for plotting
t  = tspan;

% find 1450 and 1840
indx_1450 = find(t == 1450);
indx_1840 = find(t == 1840);

% total atmospheric deposition (Mg/yr)
total_atm_dep = Matm*(k_A_tHgII + k_A_tHg0 + k_A_oHgII);

% enrichment in depisition, relative to natural
EF_atm_dep = (1/Ratm_PI(end)) * Matm;

% surface+subsurface ocean [Hg], pM
ssOceanHg = ((1e12*1e3)/(201*5.09e17))*(Mocs+Moci);

if strcmp(future,'none');
if Lplot;

    %----------------------------------------
    % Plot fast reservoirs, 1450-end year
    %----------------------------------------
    figure(81+ff);
    ff = ff+1;
    set(gca, 'FontSize',13)
    hold on
    plot (t(indx_1450:end), Matm(indx_1450:end), 'k', 'linewidth', 2.5)
    title ('Surface Hg Reservoirs')
    plot (t(indx_1450:end), Mtf(indx_1450:end),'linestyle','-.', 'linewidth',2.5,'Color',[0.7 0.7 0.7])
    plot (t(indx_1450:end), Mocs(indx_1450:end),'linestyle','--','linewidth',2.5,'color',[0.4 0.4 0.4])
    legend ('atmosphere', 'fast terrestrial', 'surface ocean', 'Location', 'NorthWest')
    xlabel('Time (years)')
    ylabel('Mg of Hg')
    xlim([t(indx_1450) t(end)])
    hold off
    
    %----------------------------------------
    % Plot intermediate and deep reservoirs,
    % 1450-2008
    %----------------------------------------
    figure(81+ff);
    ff = ff+1;
    set(gca, 'FontSize',13)
    hold on;
    plot (t(indx_1450:end), Mts(indx_1450:end), 'k', 'linewidth', 2.5)
    title ('Intermediate Hg Reservoirs')
    plot (t(indx_1450:end), Moci(indx_1450:end), 'linestyle',':', 'linewidth', 2.5,'color','k')
    plot (t(indx_1450:end), Mta(indx_1450:end), 'linestyle','--', 'linewidth', 2.5,'color',[0.5 0.5 0.5])
    plot (t(indx_1450:end), Mocd(indx_1450:end), 'linestyle','-.', 'linewidth', 2.5,'color',[0.3 0.3 0.3])
    legend('slow terrestrial', 'intermediate ocean','armored terrestrial', ...
        'deep ocean','Location', 'NorthWest')
    xlabel('Time (years)')
    ylabel('Mg of Hg')
    xlim([t(indx_1450) t(end)])
    hold off;
         
end
end

%--------------------------------------------------------------------------
% Display final mass budgets, rates, and enrichment factors
%--------------------------------------------------------------------------
  
% Display anthropogenic enrichment factors? 
L_EF = 1;  % 1 = yes {DEFAULT}
           % 0 = no

% Print output to command window
if strcmp(scenario, 'mid')
    message08 = '2008: MID-RANGE ANTHROPOGENIC EMISSION SCENARIO ';
elseif strcmp(scenario, 'high')
    message08 = '2008: HIGH ANTHROPOGENIC EMISSION SCENARIO ';
else
    message08 = '2008: LOW ANTHROPOGENIC EMISSION SCENARIO ';
end

% print model output to command window
if Ldisp;
    disp('*******************************************************************')
    disp(  message08                                                          )
    disp('*******************************************************************')
    disp(' ')
    forWeb_display_output

end

