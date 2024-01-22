function [coeffs_emis] = main_function(k_factors, Lplot, Ldisp, ...
    Lriver_FHgP, IHgD_pristine, IHgP_pristine, dt, t_SF, ...
    river_HgP_MgYr_save, river_HgD_MgYr_save)
    %==========================================================================
    % OBJECTIVE
    %   include main script in function, so that doesn't violate
    %   parallelization problems

    %--------------------------------------------------------------------------
    % Pre-anthropogenic simulation
    %--------------------------------------------------------------------------
    M = forWeb_RunPreAnthro(k_factors, Lplot, Ldisp, Lriver_FHgP, ...
        IHgD_pristine, IHgP_pristine, dt);
    % Save steady-state preindustrial reservoir sizes (Mg)
    Ratm_PI = M(1,end); % atmosphere 
    Rtf_PI  = M(2,end);  % fast terrestrial
    Rts_PI  = M(3,end);  % slow terrestrial
    Rta_PI  = M(4,end);  % armored terrestrial
    Rocs_PI = M(5,end); % surface ocean
    Roci_PI = M(6,end); % intermediate ocean
    Rocd_PI = M(7,end); % deep ocean

    R_PI = [Ratm_PI; Rtf_PI; Rts_PI; Rta_PI; Rocs_PI; Roci_PI; Rocd_PI];

    %Lplot = 1;             
    % Calculate pulse or steady emissions scenario (Mg/yr)

    % Run steady scenario
    Lpulse = 'steady';

    % figure counter
    ff = 0;
    [M_steady, ~, ~, ~] = forWeb_RunPulse(k_factors, ...
        Lplot, Ldisp, Lriver_FHgP,...
        IHgD_pristine, IHgP_pristine, R_PI, dt, ff, Lpulse, t_SF, ...
        river_HgP_MgYr_save, river_HgD_MgYr_save);

    % Run pulse scenario 
    Lpulse = 'riverpulse';
    [M_pulse, pulse_size, river_pulse,  pulse_time] = forWeb_RunPulse(k_factors,  ...
        Lplot, Ldisp, Lriver_FHgP,...
        IHgD_pristine, IHgP_pristine, R_PI, dt, ff, Lpulse, t_SF, ...
        river_HgP_MgYr_save, river_HgD_MgYr_save);

    % Analyze pulse for EAMD/EAME equations
    forWeb_AnalyzePulse
end
