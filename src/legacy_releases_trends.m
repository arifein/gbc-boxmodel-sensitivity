%  See how trajectory of releases influences trend in legacy emissions
%% 1) Analyze Streets inventory from Streets et al. (2019) all time releases paper
streets19= readmatrix('Streets2019_alltime_1500_2010_emiss_release.csv'); % air emissions
streets19_release = streets19(:,3);
streets19_emiss = streets19(:,2);
streets19_time = streets19(:,1);
%streets19_emiss(47) = streets19_emiss(47)*10;
% interpolate to annual resolution
Time    = 1510:2010;
streets19_e_i  = pchip(streets19_time, streets19_emiss , Time);  
streets19_r_i  = pchip(streets19_time, streets19_release , Time);  

% calculate additional legacy emissions from this emissions trend
time_leg = 1510:2010; % extend by 10 years for legacy emissions
leg_streets_e_19_cat = NaN(length(Time), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
leg_streets_r_19_cat = NaN(length(Time), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X

for i = 1:length(Time)
    t_p = time_leg - Time(i); % time since pulse
    t_p = t_p(t_p>0); % only select times over 0
    if (Time(i) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
        leg_streets_e_19_cat(i, i+1:end) = (0.02875 * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505 * t_p)) * streets19_e_i(i);
    else
        leg_streets_e_19_cat(i, i+1:end) = (0.0257 * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543 * t_p)) * streets19_e_i(i);
    end
    % pulse equation for releases
    leg_streets_r_19_cat(i, i+1:end) = (0.008132 * exp(-0.013 * t_p) + 0.06037 * exp(-1.369 * t_p)) * streets19_r_i(i); % 0.043% to dissolved
    %leg_streets_r_19_cat(i, i+1:end) = (0.01759 * exp(-0.01299 * t_p) + 0.1289 * exp(-1.369 * t_p)) * streets19_r_i(i); % 50% to dissolved
end

% sum these up to get total legacy emissions from air emissions
leg_streets_e_19_total = sum(leg_streets_e_19_cat,1, "omitnan");

% sum these up to get total legacy emissions from releases
leg_streets_r_19_total = sum(leg_streets_r_19_cat,1, "omitnan");

% total legacy emissions
leg_streets_19_total = leg_streets_e_19_total + leg_streets_r_19_total;
leg_streets_19_total_cat = leg_streets_e_19_cat + leg_streets_r_19_cat;

% get contributions from different time periods
leg_1510_1700 = sum(leg_streets_19_total_cat(1:191,:),1, "omitnan");
leg_1701_1800 = sum(leg_streets_19_total_cat(192:291,:),1, "omitnan");
leg_1801_1850 = sum(leg_streets_19_total_cat(292:341,:),1, "omitnan");
leg_1851_1900 = sum(leg_streets_19_total_cat(342:391,:),1, "omitnan");
leg_1901_1950 = sum(leg_streets_19_total_cat(392:441,:),1, "omitnan");
leg_1951_1990 = sum(leg_streets_19_total_cat(442:481,:),1, "omitnan");
leg_1991_2010 = sum(leg_streets_19_total_cat(482:end,:),1, "omitnan");
% stack these variables
stacked_Y = [leg_1510_1700; leg_1701_1800; leg_1801_1850; leg_1851_1900; leg_1901_1950; leg_1951_1990; leg_1991_2010];

% Make plots
figure('Position', [100 100 1000 500])
subplot(1,2,1)
plot(time_leg, leg_streets_19_total, 'linewidth',3)
hold on
plot(time_leg, leg_streets_e_19_total, 'linewidth',2)
plot(time_leg, leg_streets_r_19_total, ':', 'linewidth',2)
plot(streets19_time,streets19_emiss, '-k', 'linewidth',2)
plot(streets19_time,streets19_release, ':k', 'linewidth',2)

xlabel('Time (years)')
ylabel('Mg yr^{-1}')
legend('Legacy emissions', 'Legacy emissions from air', ...
    'Legacy emissions from releases', 'Streets19 anthro emissions',...
    'Streets19 releases', 'Location','Northwest')
axis tight
set(gca,'Fontsize',15)

subplot(1,2,2)
area(time_leg, stacked_Y')
xlabel('Time (years)')
ylabel('Mg yr^{-1}')
legend('1510-1699 legacy', '1700-1800 legacy','1801-1850 legacy','1851-1900 legacy','1901-1950 legacy', ...
    '1951-1990 legacy','1991-2010 legacy', 'Location','Northwest')
axis tight
set(gca,'Fontsize',15)

%% Keep 2010-2020 flat, increase 1970s releases by large amount, see effect on trend
% relative changes to 1970s releases
scale_min = -50; %-50% decline in 1970s releases
scale_max = 900; %up to 6 times more 1970s releases
time_leg = 2010:2020; %  years for legacy emissions
Time2 = 1510:2020; % extend by 10 years for legacy emission calculation

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at

% initialize variables
release_diff_1970 = zeros(n_runs,1); % anthropogenic release difference in 1970
leg_diff_2010_2020_1970 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_1970 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % calculate annual releases for specified scaling in 1970
    streets_release_temp = streets19_release;
    streets_release_temp(47) = streets_release_temp(47)*(100 + scale_use(i))/100; % change 1970 value
    % interpolate to annual resolution
    streets19_1970  = pchip(streets19_time, streets_release_temp , Time);  

    release_2010_2020 = repelem(streets19_1970(end), 10); % keep constant
    release_all = [streets19_1970 release_2010_2020]; % concatenate full releases timeseries
    release_diff_1970(i) = streets_release_temp(47) - streets19_release(47); % change in 1970 value

    emiss_2010_2020 = repelem(streets19_emiss(end), 10); % keep constant emissions
    emiss_all = [streets19_e_i emiss_2010_2020]; % concatenate full emissions timeseries

    leg_2010_2020_a_e = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    leg_2010_2020_a_r = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant releases at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a_e(j, :) = (0.02875 * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505 * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a_e(j, :) = (0.0257 * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543 * t_p)) * emiss_all(j);
        end
        % pulse equation for releases
        %leg_2010_2020_a_r(j, :) = (0.008132 * exp(-0.013 * t_p) + 0.06037 * exp(-1.369 * t_p)) * release_all(j); % 0.043% to dissolved
        leg_2010_2020_a_r(j, :) = (0.01759 * exp(-0.01299 * t_p) + 0.1289 * exp(-1.369 * t_p)) * release_all(j); % 50% to dissolved

    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a_e + leg_2010_2020_a_r, 1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_1970(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_1970(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(release_diff_1970, leg_diff_2010_2020_1970, 'linewidth',3)
hold on
yline(0, '--k')
xline(-1377, 'g')
xline(2938, 'g')
xlabel('Difference 1970 anthro releases compared to Streets19 (Mg yr^{-1})')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
legend('Calculation','zero line', 'P10-P90 from Streets19 for 1970')
set(gca,'Fontsize',15)
%% Keep 2010-2020 flat, test parameters: beta1
% relative changes to 1970s releases
scale_min = -2; %1/4x parameter
scale_max = 2; %4 x parameter
time_leg = 2010:2020; %  years for legacy emissions
Time2 = 1510:2020; % extend by 10 years for legacy emission calculation

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at
scale_use_fac = 2.^scale_use; % actual factor with which to scale parameter

% initialize variables
leg_diff_2010_2020_1970 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_1970 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % keep emiss and release time series as default Streets values
    release_2010_2020 = repelem(streets19_release(end), 10); % keep constant
    release_all = [streets19_r_i release_2010_2020]; % concatenate full releases timeseries

    emiss_2010_2020 = repelem(streets19_emiss(end), 10); % keep constant emissions
    emiss_all = [streets19_e_i emiss_2010_2020]; % concatenate full emissions timeseries

    leg_2010_2020_a_e = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    leg_2010_2020_a_r = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant releases at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a_e(j, :) = (0.02875 * exp(-0.01445 * scale_use_fac(i) * t_p) + 0.5268 * exp(-1.505 * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a_e(j, :) = (0.0257 * exp(-0.01678 * scale_use_fac(i) * t_p) + 0.5403 * exp(-1.543 * t_p)) * emiss_all(j);
        end
        % pulse equation for releases
        leg_2010_2020_a_r(j, :) = (0.008132 * exp(-0.013 * scale_use_fac(i) * t_p) + 0.06037 * exp(-1.369 * t_p)) * release_all(j); % 0.043% to dissolved
        %leg_2010_2020_a_r(j, :) = (0.01759 * exp(-0.01299 * t_p) + 0.1289 * exp(-1.369 * t_p)) * release_all(j); % 50% to dissolved

    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a_e + leg_2010_2020_a_r, 1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_1970(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_1970(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(scale_use_fac, leg_diff_2010_2020_1970, 'linewidth',3)
hold on
yline(0, '--k')
xlabel('Scaling factor applied to \beta_1')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
set(gca,'Fontsize',15)
%% Keep 2010-2020 flat, test parameters: beta2
% relative changes to 1970s releases
scale_min = -2; %1/4x parameter
scale_max = 2; %4 x parameter
time_leg = 2010:2020; %  years for legacy emissions
Time2 = 1510:2020; % extend by 10 years for legacy emission calculation

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at
scale_use_fac = 2.^scale_use; % actual factor with which to scale parameter

% initialize variables
leg_diff_2010_2020_1970 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_1970 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % keep emiss and release time series as default Streets values
    release_2010_2020 = repelem(streets19_release(end), 10); % keep constant
    release_all = [streets19_r_i release_2010_2020]; % concatenate full releases timeseries

    emiss_2010_2020 = repelem(streets19_emiss(end), 10); % keep constant emissions
    emiss_all = [streets19_e_i emiss_2010_2020]; % concatenate full emissions timeseries

    leg_2010_2020_a_e = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    leg_2010_2020_a_r = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant releases at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a_e(j, :) = (0.02875 * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505 * scale_use_fac(i) * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a_e(j, :) = (0.0257 * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543 * scale_use_fac(i) * t_p)) * emiss_all(j);
        end
        % pulse equation for releases
        leg_2010_2020_a_r(j, :) = (0.008132 * exp(-0.013  * t_p) + 0.06037 * exp(-1.369 * scale_use_fac(i) * t_p)) * release_all(j); % 0.043% to dissolved
        %leg_2010_2020_a_r(j, :) = (0.01759 * exp(-0.01299 * t_p) + 0.1289 * exp(-1.369 * t_p)) * release_all(j); % 50% to dissolved

    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a_e + leg_2010_2020_a_r, 1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_1970(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_1970(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(scale_use_fac, leg_diff_2010_2020_1970, 'linewidth',3)
hold on
yline(0, '--k')
xlabel('Scaling factor applied to \beta_2')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
set(gca,'Fontsize',15)
%% Keep 2010-2020 flat, test parameters: lambda1 (emiss)
% relative changes to 1970s releases
scale_min = -2; %1/4x parameter
scale_max = 2; %4 x parameter
time_leg = 2010:2020; %  years for legacy emissions
Time2 = 1510:2020; % extend by 10 years for legacy emission calculation

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at
scale_use_fac = 2.^scale_use; % actual factor with which to scale parameter

% initialize variables
leg_diff_2010_2020_1970 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_1970 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % keep emiss and release time series as default Streets values
    release_2010_2020 = repelem(streets19_release(end), 10); % keep constant
    release_all = [streets19_r_i release_2010_2020]; % concatenate full releases timeseries

    emiss_2010_2020 = repelem(streets19_emiss(end), 10); % keep constant emissions
    emiss_all = [streets19_e_i emiss_2010_2020]; % concatenate full emissions timeseries

    leg_2010_2020_a_e = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    leg_2010_2020_a_r = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant releases at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a_e(j, :) = (0.02875* scale_use_fac(i) * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505  * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a_e(j, :) = (0.0257* scale_use_fac(i) * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543  * t_p)) * emiss_all(j);
        end
        % pulse equation for releases
        leg_2010_2020_a_r(j, :) = (0.008132 * exp(-0.013  * t_p) + 0.06037 * exp(-1.369  * t_p)) * release_all(j); % 0.043% to dissolved
        %leg_2010_2020_a_r(j, :) = (0.01759 * exp(-0.01299 * t_p) + 0.1289 * exp(-1.369 * t_p)) * release_all(j); % 50% to dissolved

    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a_e + leg_2010_2020_a_r, 1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_1970(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_1970(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(scale_use_fac, leg_diff_2010_2020_1970, 'linewidth',3)
hold on
yline(0, '--k')
xlabel('Scaling factor applied to \lambda_1 (emiss)')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
set(gca,'Fontsize',15)
%% Keep 2010-2020 flat, test parameters: lambda2 (emiss)
% relative changes to 1970s releases
scale_min = -2; %1/4x parameter
scale_max = 2; %4 x parameter
time_leg = 2010:2020; %  years for legacy emissions
Time2 = 1510:2020; % extend by 10 years for legacy emission calculation

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at
scale_use_fac = 2.^scale_use; % actual factor with which to scale parameter

% initialize variables
leg_diff_2010_2020_1970 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_1970 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % keep emiss and release time series as default Streets values
    release_2010_2020 = repelem(streets19_release(end), 10); % keep constant
    release_all = [streets19_r_i release_2010_2020]; % concatenate full releases timeseries

    emiss_2010_2020 = repelem(streets19_emiss(end), 10); % keep constant emissions
    emiss_all = [streets19_e_i emiss_2010_2020]; % concatenate full emissions timeseries

    leg_2010_2020_a_e = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    leg_2010_2020_a_r = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant releases at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a_e(j, :) = (0.02875 * exp(-0.01445 * t_p) + 0.5268* scale_use_fac(i) * exp(-1.505  * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a_e(j, :) = (0.0257 * exp(-0.01678 * t_p) + 0.5403* scale_use_fac(i) * exp(-1.543  * t_p)) * emiss_all(j);
        end
        % pulse equation for releases
        leg_2010_2020_a_r(j, :) = (0.008132 * exp(-0.013  * t_p) + 0.06037 * exp(-1.369  * t_p)) * release_all(j); % 0.043% to dissolved
        %leg_2010_2020_a_r(j, :) = (0.01759 * exp(-0.01299 * t_p) + 0.1289 * exp(-1.369 * t_p)) * release_all(j); % 50% to dissolved

    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a_e + leg_2010_2020_a_r, 1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_1970(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_1970(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(scale_use_fac, leg_diff_2010_2020_1970, 'linewidth',3)
hold on
yline(0, '--k')
xlabel('Scaling factor applied to \lambda_2 (emiss)')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
set(gca,'Fontsize',15)
%% Keep 2010-2020 flat, test parameters: lambda1 (releases)
% relative changes to 1970s releases
scale_min = -2; %1/4x parameter
scale_max = 2; %4 x parameter
time_leg = 2010:2020; %  years for legacy emissions
Time2 = 1510:2020; % extend by 10 years for legacy emission calculation

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at
scale_use_fac = 2.^scale_use; % actual factor with which to scale parameter

% initialize variables
leg_diff_2010_2020_1970 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_1970 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % keep emiss and release time series as default Streets values
    release_2010_2020 = repelem(streets19_release(end), 10); % keep constant
    release_all = [streets19_r_i release_2010_2020]; % concatenate full releases timeseries

    emiss_2010_2020 = repelem(streets19_emiss(end), 10); % keep constant emissions
    emiss_all = [streets19_e_i emiss_2010_2020]; % concatenate full emissions timeseries

    leg_2010_2020_a_e = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    leg_2010_2020_a_r = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant releases at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a_e(j, :) = (0.02875 * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505  * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a_e(j, :) = (0.0257 * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543  * t_p)) * emiss_all(j);
        end
        % pulse equation for releases
        leg_2010_2020_a_r(j, :) = (0.008132* scale_use_fac(i) * exp(-0.013  * t_p) + 0.06037 * exp(-1.369  * t_p)) * release_all(j); % 0.043% to dissolved
        %leg_2010_2020_a_r(j, :) = (0.01759 * exp(-0.01299 * t_p) + 0.1289 * exp(-1.369 * t_p)) * release_all(j); % 50% to dissolved

    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a_e + leg_2010_2020_a_r, 1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_1970(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_1970(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(scale_use_fac, leg_diff_2010_2020_1970, 'linewidth',3)
hold on
yline(0, '--k')
xlabel('Scaling factor applied to \lambda_1 (releases)')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
set(gca,'Fontsize',15)
%% Keep 2010-2020 flat, test parameters: lambda2 (releases)
% relative changes to 1970s releases
scale_min = -2; %1/4x parameter
scale_max = 2; %4 x parameter
time_leg = 2010:2020; %  years for legacy emissions
Time2 = 1510:2020; % extend by 10 years for legacy emission calculation

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at
scale_use_fac = 2.^scale_use; % actual factor with which to scale parameter

% initialize variables
leg_diff_2010_2020_1970 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_1970 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % keep emiss and release time series as default Streets values
    release_2010_2020 = repelem(streets19_release(end), 10); % keep constant
    release_all = [streets19_r_i release_2010_2020]; % concatenate full releases timeseries

    emiss_2010_2020 = repelem(streets19_emiss(end), 10); % keep constant emissions
    emiss_all = [streets19_e_i emiss_2010_2020]; % concatenate full emissions timeseries

    leg_2010_2020_a_e = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    leg_2010_2020_a_r = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant releases at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a_e(j, :) = (0.02875 * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505  * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a_e(j, :) = (0.0257 * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543  * t_p)) * emiss_all(j);
        end
        % pulse equation for releases
        leg_2010_2020_a_r(j, :) = (0.008132 * exp(-0.013  * t_p) + 0.06037 * scale_use_fac(i) * exp(-1.369  * t_p)) * release_all(j); % 0.043% to dissolved
        %leg_2010_2020_a_r(j, :) = (0.01759 * exp(-0.01299 * t_p) + 0.1289 * exp(-1.369 * t_p)) * release_all(j); % 50% to dissolved

    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a_e + leg_2010_2020_a_r, 1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_1970(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_1970(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(scale_use_fac, leg_diff_2010_2020_1970, 'linewidth',3)
hold on
yline(0, '--k')
xlabel('Scaling factor applied to \lambda_2 (releases)')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
set(gca,'Fontsize',15)
%% Keep 2010-2020 flat, beta1*2 increase 1970s releases by large amount, see effect on trend
% relative changes to 1970s releases
scale_min = -50; %-50% decline in 1970s releases
scale_max = 300; %up to 6 times more 1970s releases
time_leg = 2010:2020; %  years for legacy emissions
Time2 = 1510:2020; % extend by 10 years for legacy emission calculation

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at

% initialize variables
release_diff_1970 = zeros(n_runs,1); % anthropogenic release difference in 1970
leg_diff_2010_2020_1970 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_1970 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % calculate annual releases for specified scaling in 1970
    streets_release_temp = streets19_release;
    streets_release_temp(47) = streets_release_temp(47)*(100 + scale_use(i))/100; % change 1970 value
    % interpolate to annual resolution
    streets19_1970  = pchip(streets19_time, streets_release_temp , Time);  

    release_2010_2020 = repelem(streets19_1970(end), 10); % keep constant
    release_all = [streets19_1970 release_2010_2020]; % concatenate full releases timeseries
    release_diff_1970(i) = streets_release_temp(47) - streets19_release(47); % change in 1970 value

    emiss_2010_2020 = repelem(streets19_emiss(end), 10); % keep constant emissions
    emiss_all = [streets19_e_i emiss_2010_2020]; % concatenate full emissions timeseries

    leg_2010_2020_a_e = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    leg_2010_2020_a_r = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant releases at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a_e(j, :) = (0.02875 * exp(-0.01445 *2* t_p) + 0.5268 * exp(-1.505 * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a_e(j, :) = (0.0257 * exp(-0.01678 *2* t_p) + 0.5403 * exp(-1.543 * t_p)) * emiss_all(j);
        end
        % pulse equation for releases
        leg_2010_2020_a_r(j, :) = (0.008132 *2 * exp(-0.013 * 2*t_p) + 0.06037 * exp(-1.369 * t_p)) * release_all(j); % 0.043% to dissolved
        %leg_2010_2020_a_r(j, :) = (0.01759   * exp(-0.01299 * t_p) + 0.1289 * exp(-1.369 * t_p)) * release_all(j); % 50% to dissolved

    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a_e + leg_2010_2020_a_r, 1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_1970(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_1970(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(release_diff_1970, leg_diff_2010_2020_1970, 'linewidth',3)
hold on
yline(0, '--k')
xline(-1377, 'g')
xline(2938, 'g')
xlabel('Difference 1970 anthro releases compared to Streets19 (Mg yr^{-1})')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
legend('Calculation','zero line', 'P10-P90 from Streets19 for 1970')
set(gca,'Fontsize',15)
