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
dt = 0.1;               % timestep, years

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
n_params = 40; % number of uncertainty parameters
n_samples = 1000; % number of samples to take
% initializing coefficients to store
a1 = zeros(n_samples, 1);
b1 = zeros(n_samples, 1);
a2 = zeros(n_samples, 1);
b2 = zeros(n_samples, 1);

lhs_values = lhsdesign(n_samples,n_params) * 2 -1; % create Latin Hypercube sample - between -1 and 1
f_vary = 2; % range multiplication factor between 1/f_vary and f_vary
factor_values = f_vary.^lhs_values; % calculate actual perturbation multiplication factors
tic
parfor uu = 1:n_samples
    "Uncertainty simulation " + uu
    k_factors = factor_values(uu,:); % factors to perturb parameters
    coeffs_emis = main_function(k_factors, Lplot, Ldisp, ...
        Lriver_FHgP, IHgD_pristine, IHgP_pristine, dt, t_SF, ...
        river_HgP_MgYr_save, river_HgD_MgYr_save)
    % save out results
    a1(uu) = coeffs_emis(1);
    b1(uu) = coeffs_emis(2);
    a2(uu) = coeffs_emis(3);
    b2(uu) = coeffs_emis(4);
    disp(coeffs_emis)
end
toc

save("atm_pulse_1000_sens.mat","a1","a2","b1","b2", "factor_values")