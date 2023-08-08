clear all; close all;

%==========================================================================
% VERSION 2.0
%
% OBJECTIVE
%   This module is the driver for the global biogeochemical box model of Hg
%   cycling. 
%
% CITATION FOR CODE
%   Amos, H. M., et al. (2013), Legacy impacts of all-time anthropogenic 
%     emissions on the global mercury cycle, Glob. Biogeochem. Cycle, 
%     27(2), 410-421.
%   
%   Amos, H. M., et al. (2014), Global biogeochemical implications of 
%     mercury discharges from rivers and sediment burial, Environ. Sci.
%     Technol., 48(16), 9514-9522.
%
% CREDIT
%   Co-authorship is appropriate if your paper benefits significantly from 
%   use of this model/code.
%
%   Citation is appropriate if use of this model/code has only a marginal 
%   impact on your work or if the work is a second generation application 
%   of the model/code.     
%
% NOTES
%   If you find a bug, please report it to Helen Amos 
%   (hamos@hsph.harvard.edu) with the subject "Hg box model: bug report".
%   We'll fix it, document it, and post corrected code online. 
%
%   If you'd like to submit a science update, please contact me 
%   (hamos@hsphs.harvard.edu) with the subject "Hg box model: science 
%   update". In the email, provide a quick description of the update and a 
%   copy of the journal article associated with the update. We will merge 
%   the update into the standard version of the code available online.  
%
%   We're excited to include science updates from the community, but as a 
%   policy we will ONLY include work associated with a published/accepted
%   manuscript. 
%
%   Thanks for your participation and interest!
%
% REVISION HISTORY
%   19 Jul 2012 - HMA - copied from version 5, remove the deep mineral
%                       reservoir and treat geogenic emissions as an 
%                       external forcing 
%   06 Feb 2013 - HMA - add logical for writing output to .txt file
%   30 Jan 2014 - HMA - clean up main.m driver and make RunAnthro.m into a
%                       module instead of a function, it makes the code 
%                       easier to use
%   19 Mar 2014 - HMA - add Horowitz emissions from commercial use of Hg
%   02 Sep 2014 - HMA - add option to decrease historical mining emissions
%                       (Engstrom et al., 2014)  
%   08 Sep 2014 - HMA - clean up code and comments for public release
%
% CONTACT
%   Dr. Helen M. Amos
%   Harvard Universiry
%   29 Oxford St.
%   Cambridge, MA 02138, USA
%   Email: hamos@hsph.harvard.edu
%   Phone: +1 (617)496-5348
%   Website: http://people.fas.harvard.edu/~amos/Welcome.html
%==========================================================================

%--------------------------------------------------------------------------
% Set logicals
%--------------------------------------------------------------------------

% Display plots? 
Lplot      = 0;         % 1 = yes {DEFAULT}
                        % 0 = no (faster, can be helpful for debugging) 

% Print output to command window
Ldisp      = 0;         % 1 = yes {DEFAULT}
                        % 0 = no                   
             
%--------------------------------------------------------------------------
% Set simulation time step
%--------------------------------------------------------------------------
dt = 0.2;               % timestep, years

%--------------------------------------------------------------------------
% Pick options for updated river parameterization
% Reference: Amos et al. (2014)
%--------------------------------------------------------------------------
    
% River HgD and HgP concentrations, based on observations
Lriver      = 'best';   % 'best'  = mean estimate {DEFAULT}
                        % 'low'   = mean - 1 SE
                        % 'high'  = mean + 1 SE
    
% Fraction of Hg(P) reaching open marine waters
Lriver_FHgP = 'Walsh';  % 'Walsh'   = 28% based on Zhang et al. (2014) use of
                        %             Walsh & Nitrrouer 2009 {DEFAULT}
                        % 'Chester' = 10%  based on Chester (2003)
    
% Historical scaling factors for river inputs for 
% 1970s to present. Pre-1960 scaled by releases to land and water 
% from Horowitz et al. (2014).
Lscale      = 'best';   % 'best' = best estimate {DEFAULT}
                        % 'low'  = low estimate
                        % 'high' = high estimate
                                                  
% Discharge of Hg from rivers based on Amos et al. (2014)
[IHgD_pristine, IHgP_pristine, t_SF, ... 
    river_HgP_MgYr_save, river_HgD_MgYr_save] = forWeb_riverDischarge(Lriver, Lscale, Ldisp);
%--------------------------------------------------------------------------
% Set factors for uncertainty analysis
%--------------------------------------------------------------------------
k_factors = 1;

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

%--------------------------------------------------------------------------
% Anthropogenic simulation
%--------------------------------------------------------------------------

% Select a historial emission inventory
emissInven = 'Pulse'; % 'Horowitz'  = Horowitz et al. (2014), which is 
                         %               basically Streets et al. (2011) 
                         %               plus additional emissions from 
                         %               commercial Hg use. {DEFAULT}
                         % 'Streets'   = Streets et al. (2011)
                         % 'Engstrom'  = this is the Streets et al. (2011) 
                         %               inventory with historical mining
                         %               emissions cut by 50%, as proposed 
                         %               by Engstrom et al. (2014)
                         
% Run to 2050?       
future      = 'none';    % 'none' = stop at 2008 {DEFAULT}
                         %
                         % Future scenarios from Amos et al. (2013, GBC):
                         % 'A1B'      = business-as-usual from Streets et al., 2009
                         % 'constant' = constant emisions, effectively B1 from Streets et al., 2009
                         % 'controls' = 2050 emissions are 50% of 2008 emissions
                         % 'zero'     = zero future primary anthropgenic emissions
                         
 
% Streets et al. (2011) all-time anthropogenic emissions                    
scenario    = 'mid';     % 'mid'  = Streets' central estimate {DEFAULT}
                         % 'low'  = Streets' lower 80% confidence interval
                         % 'high' = Streets' upper 80% confidence interval
%Lplot = 1;             
if (strcmp(emissInven, 'Pulse'))
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
    
    % save out results
    a1 = coeffs_emis(1);
    b1 = coeffs_emis(2);
    a2 = coeffs_emis(3);
    b2 = coeffs_emis(4);
else
    % Run with all-time anthropogenic emissions                     
    M = forWeb_RunAnthro(k_factors, Lplot, Ldisp, Lriver_FHgP, ...
        IHgD_pristine, IHgP_pristine, R_PI, dt, t_SF, ...
        river_HgP_MgYr_save, river_HgD_MgYr_save, ... 
        emissInven, future, scenario);
end
