% Analyze the results of the parameter perturbation of the Amos box model
% load('atm_pulse_1000_sens.mat')
factor_names = {'k\_A\_oHgII','k\_A\_tHgII','k\_A\_tHg0','k\_A\_oHg0', 'fdep\_tf', ...
    'fdep\_ts', 'fdep\_ta','k\_Oc\_ev','k\_Oc\_sp1','k\_Oc\_vsi', 'k\_Oc\_sp2', ...
    'k\_Oc\_vis', 'k\_Oc\_vid','k\_Oc\_sp3','k\_Oc\_vdi','k\_Te\_rf','k\_Te\_p', ...
    'k\_T\_exfs', 'k\_T\_exfa','k\_Te\_BBf\_1', 'k\_Te\_BBf\_2', 'k\_Te\_rs', ...
    'k\_T\_exsf', 'k\_T\_exsa', 'k\_Te\_BBs\_1', 'k\_Te\_BBs\_2', 'k\_Te\_ra', ...
    'k\_T\_exaf', 'k\_T\_exam', 'k\_Te\_BBa\_1', 'k\_Te\_BBa\_2', 'k\_T\_riv\_f', ...
    'k\_T\_riv\_s', 'k\_T\_riv\_a', 'k\_O\_riv\_f', 'k\_O\_riv\_s', 'k\_O\_riv\_a', ... 
    'E\_geo', 'f\_HgPexport', 'f\_diss\_r'};

%%
hist(1./b2)
%% Calculate fraction in short vs. long lifetime
short_a = (1./b1.*a1); % area under short exponential
long_a = (1./b2.*a2); % area under short exponential
short_f = short_a./(short_a + long_a); % fraction of short exponential
long_f = long_a./(short_a + long_a); % fraction of long exponential

%% short lifetime, b1
 figure('Position', [100 100 1800 1100])
 tiledlayout(8,5,'Padding', 'none', 'TileSpacing', 'compact');
for i=1:40
    nexttile
    plot(factor_values(:,i), 1./b1,'.')
    hold on
    % for moving mean, order values
    [of,ind_sort] = sort(factor_values(:,i));
    b1_sort = b1(ind_sort); % sort by x values
    b1_movemean_i = movmean(1./b1_sort,25);
    plot(of, b1_movemean_i,'-', 'LineWidth',4)
    title(factor_names{i})
    [RHO,PVAL] = corr(factor_values(:,i),1./b1,'Type','Spearman');
    text(1.7, 1.7, strcat('\rho=', num2str(round(RHO,2))), 'fontsize',13)
    if abs(RHO)>0.3
        text(0.55, 1.8, 'X', 'fontsize',15,'fontweight','bold', 'color', '#2ca25f')
    elseif abs(RHO)>0.1
        text(0.55, 1.8, 'X', 'fontsize',15, 'color', '#fdbb84')
    end
    ylim([0,2])
    xlabel('Factor scale')
    ylabel('b1 (yr)')
end
%% long lifetime, b2
 figure('Position', [100 100 1600 900])
 tiledlayout(8,5,'Padding', 'none', 'TileSpacing', 'compact');
for i=1:40
    nexttile
    plot(factor_values(:,i), 1./b2,'.')
    hold on
    % for moving mean, order values
    [of,ind_sort] = sort(factor_values(:,i));
    b2_sort = b2(ind_sort); % sort by x values
    b2_movemean_i = movmean(1./b2_sort,25);
    plot(of, b2_movemean_i,'-', 'LineWidth',4)
    title(factor_names{i})
    [RHO,PVAL] = corr(factor_values(:,i),1./b2,'Type','Spearman');
    text(1.7, 140, strcat('\rho=', num2str(round(RHO,2)))) 
end