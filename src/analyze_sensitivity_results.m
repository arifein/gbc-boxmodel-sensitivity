% Analyze the results of the parameter perturbation of the Amos box model
medium = 'atm';
load(strcat(medium,'_pulse_1000_sens.mat'))

factor_names = {'k\_A\_oHgII','k\_A\_tHgII','k\_A\_tHg0','k\_A\_oHg0', 'fdep\_tf', ...
    'fdep\_ts', 'fdep\_ta','k\_Oc\_ev','k\_Oc\_sp1','k\_Oc\_vsi', 'k\_Oc\_sp2', ...
    'k\_Oc\_vis', 'k\_Oc\_vid','k\_Oc\_sp3','k\_Oc\_vdi','k\_Te\_rf','k\_Te\_p', ...
    'k\_T\_exfs', 'k\_T\_exfa','k\_Te\_BBf\_1', 'k\_Te\_BBf\_2', 'k\_Te\_rs', ...
    'k\_T\_exsf', 'k\_T\_exsa', 'k\_Te\_BBs\_1', 'k\_Te\_BBs\_2', 'k\_Te\_ra', ...
    'k\_T\_exaf', 'k\_T\_exam', 'k\_Te\_BBa\_1', 'k\_Te\_BBa\_2', 'k\_T\_riv\_f', ...
    'k\_T\_riv\_s', 'k\_T\_riv\_a', 'k\_O\_riv\_f', 'k\_O\_riv\_s', 'k\_O\_riv\_a', ... 
    'E\_geo', 'f\_HgPexport', 'f\_diss\_r'};

%% Calculate fraction in short vs. long lifetime
short_a = (1./b1.*a1); % area under short exponential
long_a = (1./b2.*a2); % area under short exponential
short_f = short_a./(short_a + long_a); % fraction of short exponential
long_f = long_a./(short_a + long_a); % fraction of long exponential
total_int = short_a + long_a; % total integral, represents Hg emitted
%% short lifetime, b1
f = figure('Position', [100 100 1800 1100]);
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
saveas(f,strcat('Figures/sens_',medium,'_b1.png'))
%% long lifetime, b2
f = figure('Position', [100 100 1800 1100]);
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
    text(1.7, 140, strcat('\rho=', num2str(round(RHO,2))), 'fontsize',13) 
    if abs(RHO)>0.3
        text(0.55, 140, 'X', 'fontsize',15,'fontweight','bold', 'color', '#2ca25f')
    elseif abs(RHO)>0.1
        text(0.55, 140, 'X', 'fontsize',15, 'color', '#fdbb84')
    end
    xlabel('Factor scale')
    ylabel('b2 (yr)')
end
saveas(f,strcat('Figures/sens_',medium,'_b2.png'))
%% short coeff, a1
f = figure('Position', [100 100 1800 1100]);
tiledlayout(8,5,'Padding', 'none', 'TileSpacing', 'compact');
for i=1:40
    nexttile
    plot(factor_values(:,i), a1,'.')
    hold on
    % for moving mean, order values
    [of,ind_sort] = sort(factor_values(:,i));
    a1_sort = a1(ind_sort); % sort by x values
    a1_movemean_i = movmean(a1_sort,25);
    plot(of, a1_movemean_i,'-', 'LineWidth',4)
    title(factor_names{i})
    [RHO,PVAL] = corr(factor_values(:,i),a1,'Type','Spearman');
    text(1.7, 1.2, strcat('\rho=', num2str(round(RHO,2))), 'fontsize',13) 
    if abs(RHO)>0.3
        text(0.55, 1.2, 'X', 'fontsize',15,'fontweight','bold', 'color', '#2ca25f')
    elseif abs(RHO)>0.1
        text(0.55, 1.2, 'X', 'fontsize',15, 'color', '#fdbb84')
    end
    xlabel('Factor scale')
    ylabel('a1')
    ylim([0,1.3])
end
saveas(f,strcat('Figures/sens_',medium,'_a1.png'))

%% long coeff, a2
f = figure('Position', [100 100 1800 1100]);
tiledlayout(8,5,'Padding', 'none', 'TileSpacing', 'compact');
for i=1:40
    nexttile
    plot(factor_values(:,i), a2,'.')
    hold on
    % for moving mean, order values
    [of,ind_sort] = sort(factor_values(:,i));
    a2_sort = a2(ind_sort); % sort by x values
    a2_movemean_i = movmean(a2_sort,25);
    plot(of, a2_movemean_i,'-', 'LineWidth',4)
    title(factor_names{i})
    [RHO,PVAL] = corr(factor_values(:,i),a2,'Type','Spearman');
    text(1.7, 0.18, strcat('\rho=', num2str(round(RHO,2))), 'fontsize',13) 
    if abs(RHO)>0.3
        text(0.55, 0.18, 'X', 'fontsize',15,'fontweight','bold', 'color', '#2ca25f')
    elseif abs(RHO)>0.1
        text(0.55, 0.18, 'X', 'fontsize',15, 'color', '#fdbb84')
    end
    xlabel('Factor scale')
    ylabel('a2')
    ylim([0,0.2])
end
saveas(f,strcat('Figures/sens_',medium,'_a2.png'))
%% Fraction in short lifetime
f = figure('Position', [100 100 1800 1100]);
tiledlayout(8,5,'Padding', 'none', 'TileSpacing', 'compact');
for i=1:40
    nexttile
    plot(factor_values(:,i), short_f,'.')
    hold on
    % for moving mean, order values
    [of,ind_sort] = sort(factor_values(:,i));
    short_f_sort = short_f(ind_sort); % sort by x values
    short_f_movemean_i = movmean(short_f_sort,25);
    plot(of, short_f_movemean_i,'-', 'LineWidth',4)
    title(factor_names{i})
    [RHO,PVAL] = corr(factor_values(:,i),short_f,'Type','Spearman');
    text(1.7, 0.41, strcat('\rho=', num2str(round(RHO,2))), 'fontsize',13) 
    if abs(RHO)>0.3
        text(0.55, 0.41, 'X', 'fontsize',15,'fontweight','bold', 'color', '#2ca25f')
    elseif abs(RHO)>0.1
        text(0.55, 0.41, 'X', 'fontsize',15, 'color', '#fdbb84')
    end
    xlabel('Factor scale')
    ylabel('f_{short}')
    ylim([0,0.45])
end
saveas(f,strcat('Figures/sens_',medium,'_short_f.png'))
%% Total integral
f = figure('Position', [100 100 1800 1100]);
tiledlayout(8,5,'Padding', 'none', 'TileSpacing', 'compact');
for i=1:40
    nexttile
    plot(factor_values(:,i), total_int,'.')
    hold on
    % for moving mean, order values
    [of,ind_sort] = sort(factor_values(:,i));
    total_int_sort = total_int(ind_sort); % sort by x values
    total_int_movemean_i = movmean(total_int_sort,25);
    plot(of, total_int_movemean_i,'-', 'LineWidth',4)
    title(factor_names{i})
    [RHO,PVAL] = corr(factor_values(:,i),total_int,'Type','Spearman');
    text(1.7, 4.5, strcat('\rho=', num2str(round(RHO,2))), 'fontsize',13) 
    if abs(RHO)>0.3
        text(0.55, 4.3, 'X', 'fontsize',15,'fontweight','bold', 'color', '#2ca25f')
    elseif abs(RHO)>0.1
        text(0.55, 4.3, 'X', 'fontsize',15, 'color', '#fdbb84')
    end
    xlabel('Factor scale')
    ylabel('total_{area}')
    ylim([0,5])
end
saveas(f,strcat('Figures/sens_',medium,'_total_int.png'))

%% Make histograms
f = figure;
histfit(1./b1 * 12,20,'kernel')
med_b1 = prctile(1./b1 * 12, 50);
P5_b1 = prctile(1./b1 * 12, 5);
P95_b1 = prctile(1./b1 * 12, 95);
ylabel('Count')
xlabel('b1 (month)')
text(15,120, strcat("median: ", num2str(round(med_b1,1)), " months"), 'fontsize',16)
text(15,100, strcat("P5: ", num2str(round(P5_b1,1)), " months"), 'fontsize',16)
text(15,80, strcat("P95: ", num2str(round(P95_b1,1)), " months"), 'fontsize',16)
set(gca,'Fontsize',16)
saveas(f,strcat('Figures/hist_',medium,'_b1.png'))
%%
f = figure;
histfit(1./b2,20,'kernel')
med_b2 = prctile(1./b2 , 50);
P5_b2 = prctile(1./b2 , 5);
P95_b2 = prctile(1./b2 , 95);
ylabel('Count')
xlabel('b2 (yr)')
text(90,100, strcat("median: ", num2str(round(med_b2,1)), " years"), 'fontsize',16)
text(90,90, strcat("P5: ", num2str(round(P5_b2,1)), " years"), 'fontsize',16)
text(90,80, strcat("P95: ", num2str(round(P95_b2,1)), " years"), 'fontsize',16)
set(gca,'Fontsize',16)
saveas(f,strcat('Figures/hist_',medium,'_b2.png'))
%%
f = figure;
histfit(short_f,20,'kernel')
med_short_f = prctile(short_f , 50);
P5_short_f = prctile(short_f , 5);
P95_short_f = prctile(short_f , 95);
ylabel('Count')
xlabel('f_{short} ')
text(0.6,150, strcat("median: ", num2str(round(med_short_f,2))), 'fontsize',16)
text(0.6,130, strcat("P5: ", num2str(round(P5_short_f,2))), 'fontsize',16)
text(0.6,110, strcat("P95: ", num2str(round(P95_short_f,2))), 'fontsize',16)
set(gca,'Fontsize',16)
saveas(f,strcat('Figures/hist_',medium,'_short_f.png'))
%%
f = figure;
histfit(total_int(total_int<100),20,'kernel')
med_total_int = prctile(total_int , 50);
P5_total_int = prctile(total_int , 5);
P95_total_int = prctile(total_int , 95);
ylabel('Count')
xlabel('total area ')
text(4,150, strcat("median: ", num2str(round(med_total_int,2))), 'fontsize',16)
text(4,130, strcat("P5: ", num2str(round(P5_total_int,2))), 'fontsize',16)
text(4,110, strcat("P95: ", num2str(round(P95_total_int,2))), 'fontsize',16)
set(gca,'Fontsize',16)
saveas(f,strcat('Figures/hist_',medium,'_total_int.png'))


