%==========================================================================
% OBJECTIVE
%   Process pulse or steady anthropogeic Hg emissions.
%
% REFERENCES
%   Amos, H. M., et al. (2013), Legacy impacts of all-time anthropogenic 
%     emissions on the global mercury cycle, Glob. Biogeochem. Cycle, 
%     27(2), 410-421.
%   Engstrom, D. R., et al. (2014), Atmospheric hg emissions from 
%     preindustrial gold and silver extraction in the americas: A 
%     reevaluation from lake-sediment archives, Environ. Sci. Technol., 
%     48(12), 6533-6543.
%   Horowitz, H. M., et al. (2014), Historical mercury releases from 
%     commercial products: Global environmental implications, Environ. Sci.
%     Technol., 48(17), 10242-10250.
%   Nriagu, J. O. (1993), Legacy of mercury pollution, Nature, 363(6430), 
%     589-589.
%   Streets, D. G., et al. (2009), Projections of global mercury emissions 
%     in 2050, Environ. Sci. Technol., 43(8), 2983-2988.
%   Streets, D. G., et al. (2011), All-time releases of mercury to the 
%     atmosphere from human activities, Environ. Sci. Technol., 45(24), 
%     10485-10491.   
%
% REVISION HISTORY
%   02 Nov 2011 - HMA - now interpolate with pchip() instead of spline()
%   08 Nov 2011 - HMA - use built-in MATLAB functions to optimize the code
%   14 Jan 2012 - hma - extrend pre-1850 anthropogenic emissions to 2000 BC
%   17 Jan 2012 - hma - add time-dependent volcanic eruptions
%   27 Apr 2012 - hma - add David Streets' new 2050 SRES scenarios
%   03 May 2012 - hma - extrapolate future scenarios to 2100
%   30 Jul 2012 - hma - clean up code/comments
%   30 Jan 2014 - HMA - update filepath to Streets anthro emission files
%   02 May 2014 - HMA - create figure of primary emissions from 1500 to
%                       2050 for private defense slides
%   02 Sep 2014 - HMA - add option to decrease historical mining emissions
%                       (Engstrom et al., 2014)
%   08 Sep 2014 - HMA - clean up code and comments for public release
%   28 Mar 2023 - AF - edit code for pulse emissions setup
%==========================================================================

%%
%--------------------------------------------------------------------------
% Emission inventory from Streets et al. (2011) - run as historical
%--------------------------------------------------------------------------

% read in emissions from .txt file
load('AnthroEmissAllTime_20120112.txt')

% parse data
Syear   = AnthroEmissAllTime_20120112(:,1); % 2000 BC to 2008 AD, decadal
Streets = AnthroEmissAllTime_20120112(:,2); % Mg/yr

%%

% interpolate to annual resolution
Time    = (Syear(1):Syear(end));
Anthro  = pchip(Syear, Streets , Time);  % mid


%%
%--------------------------------------------------------------------------
% Future emission scenarios - steady emissions, or pulse

clear n m; % for safety's sake
y_end = 2110;
FTime = [Time, 2009:1:y_end]; % future time, annual resoltuion
n = length(Anthro)+1;
m = y_end - 2009; % treaty comes into force in 2016
Anthro_steady = zeros(1,n+m);
Anthro_pulse = zeros(1, n+m);
Anthro_steady(1:n-1) = Anthro;
Anthro_pulse(1:n-1) = Anthro;
% initialize future scenarios with 2009
for i = n:n + m
    Anthro_steady(i)  = Anthro(end);   % unchanging
    Anthro_pulse(i)  = Anthro(end);   % unchanging
end
% generate pulse of X Mg/yr in specific year
pulse_idx = n+1; % index of time where want to input pulse - 2010
% pulse_idx = 3501; % index of time where want to input pulse - 1500
% pulse_idx = 3751; % index of time where want to input pulse - 1750
% pulse_idx = 3901; % index of time where want to input pulse - 1900
% pulse_idx = 3971; % index of time where want to input pulse - 1970
% pulse_idx = 3996; % index of time where want to input pulse - 1995
%pulse_idx = 3849; % index of time where want to input pulse - 1848

pulse_time = FTime(pulse_idx); % what year is this actually
% set pulse to 0 if want only a river pulse
pulse_size = 0; % how big is pulse in Mg/yr, results are independent
Anthro_pulse(pulse_idx) = Anthro_pulse(pulse_idx) + pulse_size; % implement pulse

    
%% 
%--------------------------------------------------------------------------
% Emissions from intentional use of Hg in commericial products and
% processes
%
% Reference: Horowitz et al. (2014)
%--------------------------------------------------------------------------

% Use the Horowitz et al. (2014) inventory, this automatically sets
% a logical that tells the code to add commericial Hg emissions to the
% Streets et al. (2011) inventory
Lprod = 1;
forWeb_product_emissions
% fill in with additional steady values so that the same size as Streets
Ep_lf = [Ep_lf repelem(Ep_lf(end), m+1)];
Ep_atm = [Ep_atm repelem(Ep_atm(end), m+1)];
E_Streets_waste = [E_Streets_waste repelem(E_Streets_waste(end), m+1)];
E_Streets_overlap = [E_Streets_overlap repelem(E_Streets_overlap(end), m+1)];

% Add on product emissions to Streets emissions
Anthro_pulse = Anthro_pulse - E_Streets_waste + Ep_atm;
Anthro_steady = Anthro_steady - E_Streets_waste + Ep_atm;

%%
%--------------------------------------------------------------------------
% PLOTS
%--------------------------------------------------------------------------

if Lplot;
    % primary anthropogenic emissions, including future scenarios
    hfig = figure(80);
    set(hfig,'units','normalized','Position',[0.1 0.4 0.5 .7])
    set(gcf,'Color',[1 1 1])
    set(gca,'FontSize',18)
    hold on;
    plot(FTime,Anthro_steady  ,'b','LineWidth',3)
    plot(FTime,Anthro_pulse,'g','LineWidth',3)
    xlabel('Year (AD)')
    ylabel('(Mg a^{-1}) ')
    xlim([1500 y_end])
    title('Primary Anthropogenic Emissions for Pulse Case')
    hold off;
    legend('Steady','Pulse','Location','NorthWest')
    
end
