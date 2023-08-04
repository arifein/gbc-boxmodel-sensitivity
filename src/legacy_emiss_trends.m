% Load different emissions inventories (Streets, EDGAR, Horowitz, etc...)
% and see how trajectory influences trend in legacy emissions
%% 1) Analyze Streets inventory from Streets et al. (2019) all time releases paper
streets19= readmatrix('Streets2019_alltime_1500_2010_emiss.csv'); % air emissions
streets19_emiss = streets19(:,2);
streets19_time = streets19(:,1);
%streets19_emiss(47) = streets19_emiss(47)*10;
% interpolate to annual resolution
Time    = 1510:2010;
streets19_i  = pchip(streets19_time, streets19_emiss , Time);  

% calculate additional legacy emissions from this emissions trend
time_leg = 1510:2010; % extend by 10 years for legacy emissions
leg_streets_19_cat = NaN(length(Time), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X

for i = 1:length(Time)
    t_p = time_leg - Time(i); % time since pulse
    t_p = t_p(t_p>0); % only select times over 0
    if (Time(i) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
        leg_streets_19_cat(i, i+1:end) = (0.02875 * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505 * t_p)) * streets19_i(i);
    else
        leg_streets_19_cat(i, i+1:end) = (0.0257 * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543 * t_p)) * streets19_i(i);
    end
end

% sum these up to get total legacy emissions
leg_streets_19_total = sum(leg_streets_19_cat,1, "omitnan");
% get contributions from different time periods
leg_1510_1699 = sum(leg_streets_19_cat(1:190,:),1, "omitnan");
leg_1700_1849 = sum(leg_streets_19_cat(191:340,:),1, "omitnan");
leg_1850_1899 = sum(leg_streets_19_cat(341:390,:),1, "omitnan");
leg_1900_1949 = sum(leg_streets_19_cat(391:440,:),1, "omitnan");
leg_1950_1979 = sum(leg_streets_19_cat(441:470,:),1, "omitnan");
leg_1980_1999 = sum(leg_streets_19_cat(471:490,:),1, "omitnan");
leg_2000_2010 = sum(leg_streets_19_cat(491:end,:),1, "omitnan");
% stack these variables
stacked_Y = [leg_1510_1699; leg_1700_1849; leg_1850_1899; leg_1900_1949; leg_1950_1979; leg_1980_1999; leg_2000_2010];

% Make plots
figure('Position', [100 100 1000 500])
subplot(1,2,1)
plot(time_leg, leg_streets_19_total, 'linewidth',3)
hold on
plot(streets19_time,streets19_emiss, 'linewidth',3)
xlabel('Time (years)')
ylabel('Mg yr^{-1}')
legend('Legacy emissions', 'Streets19 anthro emissions','Location','Northwest')
axis tight
set(gca,'Fontsize',15)

subplot(1,2,2)
area(time_leg, stacked_Y')
xlabel('Time (years)')
ylabel('Mg yr^{-1}')
legend('1510-1699 legacy', '1700-1849 legacy','1850-1899 legacy','1900-1949 legacy', ...
    '1950-1979 legacy','1980-1999 legacy', '2000-2009 legacy', 'Location','Northwest')
axis tight
set(gca,'Fontsize',15)

%% Try to vary 2010–2020 trend, what will be effect on legacy emissions?
% vary legacy trend for 2010-2020 from -90% per decade to 200% per decade
% calculate additional legacy emissions from this emissions trend
time_leg = 2010:2020; %  years for legacy emissions
Time2 = 1510:2020; % extend by 10 years for legacy emission calculation

% relative changes in emissions
trend_min = -90; %-90% decline over 2010-2020
trend_max = 100; %doubling over 2010-2020

n_runs = 50; % number of runs to choose
trend_use = linspace(trend_min, trend_max, n_runs); % values to sample at

% initialize variables
emiss_diff_2010_2020 = zeros(n_runs,1); % anthropogenic emission difference between 2010 and 2020
leg_diff_2010_2020 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(trend_use)
    % calculate annual emissions for specified trend
    emiss_2010_2020 = linspace(streets19_i(end), streets19_i(end) * (100 + trend_use(i))/100, 11);
    emiss_all = [streets19_i emiss_2010_2020(2:end)]; % concatenate full emiss timeseries
    emiss_diff_2010_2020(i) = emiss_2010_2020(end) - emiss_2010_2020(1);
    
    leg_2010_2020_a = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a(j, :) = (0.02875 * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505 * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a(j, :) = (0.0257 * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543 * t_p)) * emiss_all(j);
        end
    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a,1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(emiss_diff_2010_2020, leg_diff_2010_2020, 'linewidth',3)
hold on
yline(0, '--k')
xlabel('Difference 2020 minus 2010 anthro emiss (Mg yr^{-1})')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
set(gca,'Fontsize',15)

% subplot(1,2,2)
% plot(emiss_diff_2010_2020, leg_trend_2010_2020, 'linewidth',3)
% hold on
% yline(0, '--k')
% xlabel('Difference 2020 minus 2010 anthro emiss (Mg yr^{-1})')
% ylabel('Trend 2010–2020 legacy emiss (Mg yr^{-2})')
% axis tight
% set(gca,'Fontsize',15)

%% Keep 2010-2020 flat, increase 1970s emissions by large amount, see effect on trend
% relative changes to 1970s emissions
scale_min = -50; %-50% decline in 1970s emiss
scale_max = 500; %up to 6 times more 1970s emiss

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at

% initialize variables
emiss_diff_1970 = zeros(n_runs,1); % anthropogenic emission difference in 1970
leg_diff_2010_2020_1970 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_1970 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % calculate annual emissions for specified scaling in 1970
    streets_emiss_temp = streets19_emiss;
    streets_emiss_temp(47) = streets_emiss_temp(47)*(100 + scale_use(i))/100; % change 1970 value
    % interpolate to annual resolution
    streets19_1970  = pchip(streets19_time, streets_emiss_temp , Time);  

    emiss_2010_2020 = repelem(streets19_1970(end), 10); % keep constant
    emiss_all = [streets19_1970 emiss_2010_2020]; % concatenate full emiss timeseries
    emiss_diff_1970(i) = streets_emiss_temp(47) - streets19_emiss(47); % change in 1970 value
    
    leg_2010_2020_a = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a(j, :) = (0.02875 * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505 * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a(j, :) = (0.0257 * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543 * t_p)) * emiss_all(j);
        end
    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a,1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_1970(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_1970(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(emiss_diff_1970, leg_diff_2010_2020_1970, 'linewidth',3)
hold on
yline(0, '--k')
xlabel('Difference 1970 anthro emissions compared to Streets19 (Mg yr^{-1})')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
set(gca,'Fontsize',15)
%% Keep 2010-2020 flat, increase 2000s emissions by large amount, see effect on trend
% relative changes to 2000s emissions
scale_min = -50; %-50% decline in 2000s emiss
scale_max = 500; %up to 6 times more 2000s emiss

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at

% initialize variables
emiss_diff_2000 = zeros(n_runs,1); % anthropogenic emission difference in 2000
leg_diff_2010_2020_2000 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_2000 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % calculate annual emissions for specified scaling in 2000
    streets_emiss_temp = streets19_emiss;
    streets_emiss_temp(50) = streets_emiss_temp(50)*(100 + scale_use(i))/100; % change 2000 value
    % interpolate to annual resolution
    streets19_2000  = pchip(streets19_time, streets_emiss_temp , Time);  

    emiss_2010_2020 = repelem(streets19_2000(end), 10); % keep constant
    emiss_all = [streets19_2000 emiss_2010_2020]; % concatenate full emiss timeseries
    emiss_diff_2000(i) = streets_emiss_temp(50) - streets19_emiss(50); % change in 2000 value
    
    leg_2010_2020_a = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a(j, :) = (0.02875 * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505 * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a(j, :) = (0.0257 * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543 * t_p)) * emiss_all(j);
        end
    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a,1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_2000(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_2000(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(emiss_diff_2000, leg_diff_2010_2020_2000, 'linewidth',3)
hold on
yline(0, '--k')
xlabel('Difference 2000 anthro emissions compared to Streets19 (Mg yr^{-1})')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
set(gca,'Fontsize',15)
%% Keep 2010-2020 flat, increase 1990s emissions by large amount, see effect on trend
% relative changes to 1990s emissions
scale_min = -50; %-50% decline in 1990s emiss
scale_max = 500; %up to 6 times more 1990s emiss

n_runs = 50; % number of runs to choose
scale_use = linspace(scale_min, scale_max, n_runs); % values to sample at

% initialize variables
emiss_diff_1990 = zeros(n_runs,1); % anthropogenic emission difference in 1990
leg_diff_2010_2020_1990 = zeros(n_runs,1); % legacy emission difference between 2010 and 2020
leg_trend_2010_2020_1990 = zeros(n_runs,1); % legacy emission trend between 2010 and 2020

for i = 1:length(scale_use)
    % calculate annual emissions for specified scaling in 1990
    streets_emiss_temp = streets19_emiss;
    streets_emiss_temp(49) = streets_emiss_temp(49)*(100 + scale_use(i))/100; % change 1990 value
    % interpolate to annual resolution
    streets19_1990  = pchip(streets19_time, streets_emiss_temp , Time);  

    emiss_2010_2020 = repelem(streets19_1990(end), 10); % keep constant
    emiss_all = [streets19_1990 emiss_2010_2020]; % concatenate full emiss timeseries
    emiss_diff_1990(i) = streets_emiss_temp(49) - streets19_emiss(49); % change in 1990 value
    
    leg_2010_2020_a = NaN(length(Time2), length(time_leg)); % legacy emissions at year Y driven by ant emissions at year X
    for j = 1:length(Time2)
        t_p = time_leg - Time2(j); % time since pulse
        t_p(t_p<=0) = NaN; % for times before pulse, use NaN values
        if (Time2(j) < 1849) % use different equations for pre-ind and post-ind, from pulse expts
            leg_2010_2020_a(j, :) = (0.02875 * exp(-0.01445 * t_p) + 0.5268 * exp(-1.505 * t_p)) * emiss_all(j);
        else
            leg_2010_2020_a(j, :) = (0.0257 * exp(-0.01678 * t_p) + 0.5403 * exp(-1.543 * t_p)) * emiss_all(j);
        end
    end
    
    % calculate total legacy emissions for 2010–2020 time period
    leg_2010_2020_total = sum(leg_2010_2020_a,1, "omitnan");
    
    % calculate difference in legacy emissions for 2020 minus 2010    
    leg_diff_2010_2020_1990(i) = leg_2010_2020_total(end) - leg_2010_2020_total(1);
    
    % calculate trend
    p = polyfit(time_leg,leg_2010_2020_total,1); % slope and intercept for linear fit
    leg_trend_2010_2020_1990(i) = p(1);
end

% Make plots
figure('Position', [100 100 1000 500])
%subplot(1,2,1)
plot(emiss_diff_1990, leg_diff_2010_2020_1990, 'linewidth',3)
hold on
yline(0, '--k')
xlabel('Difference 1990 anthro emissions compared to Streets19 (Mg yr^{-1})')
ylabel('Difference 2020 minus 2010 legacy emiss (Mg yr^{-1})')
axis tight
set(gca,'Fontsize',15)

%% Statistical relationship - how to relate previous emissions amount + current trend to legacy emissions trend