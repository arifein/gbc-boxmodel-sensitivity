function M = forWeb_RunPreAnthro(k_factors, Lplot, Ldisp, Lriver_FHgP, ...
    IHgD_pristine, IHgP_pristine, dt)
    %==========================================================================
    % OBJECTIVE
    %   This module simulates natural, pre-anthropogenic steady-state
    %   conditions, which then serve as the initial condition for the
    %   anthropogenic era simulation. 
    %
    % REVISION HISTORY
    %   05 Jun 2011 - EMS - (listed as last modification date)
    %   02 Nov 2011 - HMA - strucutural changes, move things around for
    %                       organization
    %   08 Nov 2011 - HMA - utilize built-in MATLAB functions to optimize code
    %   06 Dec 2011 - hma - specify different biomass burning emissions for the
    %                       prei-1450 period and 1450-2008
    %   07 Dec 2011 - hma - add plots to check the emissions/deposition balance
    %                       of the atmosphere
    %   07 Dec 2011 - hma - bug fix: biomass burning now constant emission
    %                       fluxes and not rate constants
    %   31 Dec 2011 - hma - biomass burning now f(t), not a constant flux
    %   18 Jul 2012 - hma - plot difference between burial and geogenic
    %                       emissions
    %   19 Jul 2012 - hma - treat geogenic emissions as an external forcing and
    %                       remove deep mineral reservoir from the system
    %   07 Aug 2012 - hma - move Fate of Hg test (Figure 5, GBC paper) to
    %                       FateOfHg.m
    %==========================================================================

    % Assemble matrix A, such that dM/dt = A*M + E
    sim_type = 1;  % 1 = pre-anthropogenic, 2 = anthropogenic era
    [A, E_geo, k_mat] = forWeb_makeA(k_factors, sim_type, Lriver_FHgP, IHgD_pristine, IHgP_pristine);

    % integration time (yrs)
    tspan = 0:dt:6e4;

    % external emissions (Mg/yr)
    EE = diag(ones(1,7),0);         % indentity matrix
    Eg = [E_geo; 0; 0; 0; 0; 0; 0]; % geogenic emissions to the atmosphere

    % dummy matrix of zeros (not necessary, but dramatically saves time)
    M = zeros(7, numel(tspan));

    % Initial reservoir sizes
    % Atmosphere
    Ratm       = 4366;               % Feinberg et al. (2022)

    % Ocean
    Rocs       = 2407;               % Zhang et al. (2014)
    Roci       = 134000;             % Sunderland and Mason (2007)
    Rocd       = 220649;             % Sunderland and Mason (2007)

    %  Terrestrial
    Rtf        = 9620;               % Leaf, fast and intermediate pools from Smith-Downey et al (2010)
    Rts        = 34900;              % Slow pool from Smith-Downey et al. (2010)
    Rta        = 193600;             % Armored pool from Smith-Downey et al. (2010)

    % Initial conditions (Mg)
    M(:,1) = [Ratm ; Rtf; Rts; Rta; Rocs; Roci; Rocd]; % normal initial conditions

    % Solve M(t) forward in time
    for j = 2:numel(tspan);

        % with external forcing from geogenic emissions
        M(:,j) = ( A*M(:,j-1) + EE*Eg )*dt + M(:,j-1);

    end

    % parse output
    Matm    = M(1,:);    % atmosphere
    Mtf     = M(2,:);    % fast soil pool
    Mts     = M(3,:);    % slow soil pool
    Mta     = M(4,:);    % armored soil pool
    Mocs    = M(5,:);    % surface ocean
    Moci    = M(6,:);    % intermediate ocean
    Mocd    = M(7,:);    % deep ocean

    % Save steady-state preindustrial reservoir sizes (Mg)
    Ratm_PI = Matm(end); % atmosphere
    Rtf_PI  = Mtf(end);  % fast terrestrial
    Rts_PI  = Mts(end);  % slow terrestrial
    Rta_PI  = Mta(end);  % armored terrestrial
    Rocs_PI = Mocs(end); % surface ocean
    Roci_PI = Moci(end); % intermediate ocean
    Rocd_PI = Mocd(end); % deep ocean


    %--------------------------------------------------------------------------
    % PLOTS
    %--------------------------------------------------------------------------

    if Lplot;

        %---------------------------------------
        % Plot reservoir masses vs. time
        %---------------------------------------
        figure(1)
        plot (tspan, Matm, 'g-', 'linewidth', 2)
        title ('Fast Hg Reservoirs (Mg)')
        hold on
        plot (tspan, Mtf, 'r-', 'linewidth', 2)
        plot (tspan, Mocs, 'b-', 'linewidth', 2)
        legend ('atmosphere', 'fast terrestrial','surface ocean', 'Location', 'NorthOutside')

        figure(2)
        plot (tspan, Mts, 'g-', 'linewidth', 2)
        title ('Intermediate Hg Reservoirs (Mg)')
        hold on
        plot (tspan, Moci, 'b-', 'linewidth', 2)
        legend ('slow terrestrial', 'intermediate ocean','Location', 'NorthOutside')
        hold off

        figure(3)
        plot (tspan, Mta, 'g-', 'linewidth', 2)
        title ('Slow Hg Reservoirs (Mg)')
        hold on
        plot (tspan, Mocd, 'b-', 'linewidth', 2)
        legend ('armored terrestrial', 'deep ocean','Location', 'NorthOutside')
        hold off

        %---------------------------------------
        % plot the difference between burial and geogenic emissions
        %---------------------------------------
        figure(4)
        set(gca,'FontSize',14)
        plot(tspan,E_geo - k_Oc_sp3*Mocd)
        xlabel('Time (years)')
        ylabel('Geogenic Emissions - Burial, Mg a^{-1} ')
        title('Steady State Check')   

    end

    %--------------------------------------------------------------------------
    % Display final mass budgets and rates
    %--------------------------------------------------------------------------

    % Display enrichment factors?
    %
    % There's always 1 (i.e., zero enrichment) for the pre-anthropogenic
    % simulation, so the default is to not display them. 
    L_EF = 0;  % 1 = yes
               % 0 = no {DEFAULT}

    if Ldisp;
        % Print output to screen
        disp('*******************************************************************')
        disp('PRE-ANTHROPOGENIC STEADY STATE VALUES '                             )
        disp('*******************************************************************')
        disp(' ')

        % print output to command window
        forWeb_display_output

    end
end

