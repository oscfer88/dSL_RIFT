%% Aplha-power analysis

clear
ft_defaults

% set params
selsensors      = 'occ';

% select files
folder = 'Q:/MEG_data/';
filelist = dir(fullfile(folder, '**', 'tfr_all3.mat'));

filelistHiLo = filelist([1 5 3 8  9 13 12 16 17 20 21 24 25 29 27 34 32 36 37 40]); %HIGH-LOW condition
filelistLoHi = filelist([7 2 6 4 11 10 15 14 19 18 23 22 28 26 31 30 35 33 39 38]); %LOW-HIGH condition

% load data
nsubj = size(filelistHiLo,1);
tfr_dataHiLo = cell(nsubj,1);
tfr_dataLoHi = cell(nsubj,1);
for i = 1:nsubj
    tfr_tempHiLo = load([filelistHiLo(i).folder '\' filelistHiLo(i).name]);
    tfr_dataHiLo{i} = tfr_tempHiLo.tfr_all;
    tfr_tempLoHi = load([filelistLoHi(i).folder '\' filelistLoHi(i).name]);
    tfr_dataLoHi{i} = tfr_tempLoHi.tfr_all;
end
tfr_dataHiLo_raw = tfr_dataHiLo; %tfr_dataHiLo = tfr_dataHiLo_raw;
tfr_dataLoHi_raw = tfr_dataLoHi; %tfr_dataLoHi = tfr_dataLoHi_raw;

% baseline correction
cfg               = [];
cfg.baseline      = [-1 -0.5];
cfg.baselinetype  = 'relative'; %NOTE: use 'relchange' in the future
for i = 1:nsubj
    tfr_dataHiLo{i} = ft_freqbaseline(cfg, tfr_dataHiLo{i});
    tfr_dataLoHi{i} = ft_freqbaseline(cfg, tfr_dataLoHi{i});
end

% combine sensors
cfg                 = [];
cfg.method          = 'sum';
for i = 1:nsubj
    tfr_dataHiLo{i} = ft_combineplanar(cfg,tfr_dataHiLo{i}); %ft_combineplanar_nocheck
    tfr_dataLoHi{i} = ft_combineplanar(cfg,tfr_dataLoHi{i});
end

% difference
tfr_data_dif = cell(nsubj,1);
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
for i = 1:nsubj
    tfr_dataDiff{i} = ft_math(cfg,tfr_dataHiLo{i},tfr_dataLoHi{i});
end

% grand average
cfg                 = [];
tfr_dataHiLo_split = ft_freqgrandaverage(cfg, tfr_dataHiLo{:});
tfr_dataLoHi_split = ft_freqgrandaverage(cfg, tfr_dataLoHi{:});
tfr_dataDiff_ga = ft_freqgrandaverage(cfg, tfr_dataDiff{:});


% topoplots
plotmarkers = 1;
cfg                 = []; %plot topography of alpha (8-12 HZ)
cfg.xlim            = [0 1.5];
cfg.ylim            = [8 12];
cfg.zlim            = [.5  3];
cfg.layout  = 'neuromag306cmb.lay';
if plotmarkers
    left_sensors  = [121 122 139 140 123 124 127 128 129 130 141 142 143 144 145 146 147 148 155 156];
    right_sensors = [169 170 187 188 185 186 191 192 189 190 173 174 179 180 177 178 175 176 153 154];
    left_sensors = left_sensors([2:2:numel(left_sensors)])/2;
    right_sensors = right_sensors([2:2:numel(right_sensors)])/2;

    cfg.highlight           = 'on';
    cfg.highlightchannel    = [left_sensors, right_sensors];
    cfg.highlightsize       = 10;
end
figure; ft_topoplotTFR(cfg, tfr_dataHiLo_split); title('High-Low');
colormap jet

figure; ft_topoplotTFR(cfg, tfr_dataLoHi_split); title('Low-High');
colormap jet

cfg.zlim            = [-.15 .15];
figure; ft_topoplotTFR(cfg, tfr_dataDiff_ga); title('Difference');
colormap jet


% average over frequency band, time window and selected sensors
tfr_dataHiLo_leftpost = cell(nsubj,1);
tfr_dataHiLo_rightpost = cell(nsubj,1);
tfr_dataLoHi_leftpost = cell(nsubj,1);
tfr_dataLoHi_rightpost = cell(nsubj,1);
if strcmp(selsensors, 'occ')
    left_sensors  = [121 122 139 140 123 124 127 128 129 130 141 142 143 144 145 146 147 148 155 156];
    right_sensors = [169 170 187 188 185 186 191 192 189 190 173 174 179 180 177 178 175 176 153 154];
else
    left_sensors  = [15 16 13 14 31 32 117 118 119 120 133 134 123 124 121 122 139 140];
    right_sensors = [97 98 99 100 81 82 183 184 181 182 167 168 185 186 187 188 169 170];
end
if strcmp(combineplanar, 'yes')
    left_sensors = left_sensors([2:2:numel(left_sensors)])/2;
    right_sensors = right_sensors([2:2:numel(right_sensors)])/2;
end
cfg                 = [];
cfg.frequency       = [8 12];   %alpha
cfg.avgoverfreq     = 'yes';
cfg.latency         = [0 1.5];  %pre-search array time window
cfg.avgovertime     = 'yes';
cfg.channel         = left_sensors;
cfg.avgoverchan     = 'yes';
for i = 1:nsubj
    tfr_dataHiLo_leftpost{i} = ft_selectdata(cfg, tfr_dataHiLo{i});
    tfr_dataLoHi_leftpost{i} = ft_selectdata(cfg, tfr_dataLoHi{i});
end
cfg.channel         = right_sensors;
for i = 1:nsubj
    tfr_dataHiLo_rightpost{i} = ft_selectdata(cfg, tfr_dataHiLo{i});
    tfr_dataLoHi_rightpost{i} = ft_selectdata(cfg, tfr_dataLoHi{i});
end

% individual average
HiLo_leftpost_sbj = nan(nsubj,1);
LoHi_leftpost_sbj = nan(nsubj,1);
HiLo_rightpost_sbj = nan(nsubj,1);
LoHi_rightpost_sbj = nan(nsubj,1);
for i = 1:nsubj
    HiLo_leftpost_sbj(i) = tfr_dataHiLo_leftpost{i}.powspctrm;
    LoHi_leftpost_sbj(i) = tfr_dataLoHi_leftpost{i}.powspctrm;
    HiLo_rightpost_sbj(i) = tfr_dataHiLo_rightpost{i}.powspctrm;
    LoHi_rightpost_sbj(i) = tfr_dataLoHi_rightpost{i}.powspctrm;
end

% grand average
cfg         = [];
tfr_dataHiLo_leftpost_ga = ft_freqgrandaverage(cfg, tfr_dataHiLo_leftpost{:});
HiLo_leftpost_ga = tfr_dataHiLo_leftpost_ga.powspctrm; clearvars tfr_dataHiLo_leftpost_ga;
tfr_dataLoHi_leftpost_ga = ft_freqgrandaverage(cfg, tfr_dataLoHi_leftpost{:});
LoHi_leftpost_ga = tfr_dataLoHi_leftpost_ga.powspctrm; clearvars tfr_dataLoHi_leftpost_ga;
tfr_dataHiLo_rightpost_ga = ft_freqgrandaverage(cfg, tfr_dataHiLo_rightpost{:});
HiLo_rightpost_ga = tfr_dataHiLo_rightpost_ga.powspctrm; clearvars tfr_dataHiLo_rightpost_ga;
tfr_dataLoHi_rightpost_ga = ft_freqgrandaverage(cfg, tfr_dataLoHi_rightpost{:});
LoHi_rightpost_ga = tfr_dataLoHi_rightpost_ga.powspctrm; clearvars tfr_dataLoHi_rightpost_ga;



% Plot for Fig. 5 of the paper [part 1]
addpath(['...\Violinplot'])
% plot alpha power by condition (hemi x probability)
f1 = figure('Units', 'centimeters', 'Position', [5,5,8,8]);
subplot(221)
vp = violinplot([HiLo_rightpost_sbj, LoHi_rightpost_sbj, LoHi_leftpost_sbj, HiLo_leftpost_sbj], [], 'ShowMean', true);
xlim([0.5,4.5]); xlabel('Distractor probability'); xticklabels({'high', 'low', 'high', 'low'});
ylim([0.4 2.1]); ylabel('Alpha power');
legend('left','right','FontSize',6); legend boxoff
ax = gca;
ax.FontSize = 6;
vp(1, 1).ViolinColor = [0.3010    0.7450    0.9330];
vp(1, 2).ViolinColor = [0.3010    0.7450    0.9330];
vp(1, 3).ViolinColor = [0.8500    0.3250    0.0980];
vp(1, 4).ViolinColor = [0.8500    0.3250    0.0980];
% plot line between each subject data point
x1 = vp(1,1).ScatterPlot.XData;
y1 = vp(1,1).ScatterPlot.YData;
x2 = vp(1,2).ScatterPlot.XData;
y2 = vp(1,2).ScatterPlot.YData;
plot([x1; x2],[y1; y2],'Color',[.8 .8 .8],'linewidth',0.5)
x3 = vp(1,3).ScatterPlot.XData;
y3 = vp(1,3).ScatterPlot.YData;
x4 = vp(1,4).ScatterPlot.XData;
y4 = vp(1,4).ScatterPlot.YData;
plot([x3; x4],[y3; y4],'Color',[.8 .8 .8],'linewidth',0.5)
% add legend
hleglines = [vp(1, 1).ViolinPlot, vp(1, 3).ViolinPlot];
legend(hleglines, {'left', 'right'},'Location','northeast')
legend('boxoff')


% average across conditions
alpha_contraHi = mean([HiLo_rightpost_sbj LoHi_leftpost_sbj],2);
alpha_contraLo = mean([HiLo_leftpost_sbj LoHi_rightpost_sbj],2);

% grand average
alpha_contraHi_ga = mean(alpha_contraHi);
alpha_contraLo_ga = mean(alpha_contraLo);


% Plot for Fig. 5 of the paper [part 2]
% plot delta alpha power (dif between probabilities)
figure(f1)
subplot(222)
alpha_diff_sbj = alpha_contraLo - alpha_contraHi;
violinplot(alpha_diff_sbj, [], 'ShowMean', true); hold on
yline(0, '--'); hold off
xlim([0.5,1.5]); xticklabels('low>high');
ylim([-.1, .05]); ylabel('Alpha power');
ax = gca;
ax.FontSize = 6;


% statistic
addpath('...\Cohens_d')
addpath('...\Bayes_factor')

[~,p,~,stats] = ttest(alpha_contraHi, alpha_contraLo) %t-test
computeCohen_d(alpha_contraHi, alpha_contraLo, 'paired')
[bf10,p] = bf.ttest(alpha_contraHi,alpha_contraLo)


%% Alpha-power time course

% select sensors
ptime = tfr_dataHiLo{1}.time;

tfr_dataHiLo_left = cell(nsubj,1);
tfr_dataLoHi_left = cell(nsubj,1);
tfr_dataHiLo_right = cell(nsubj,1);
tfr_dataLoHi_right = cell(nsubj,1);
HiLo_left_sbj = nan(nsubj,size(ptime,2));
LoHi_left_sbj = nan(nsubj,size(ptime,2));
HiLo_right_sbj = nan(nsubj,size(ptime,2));
LoHi_right_sbj = nan(nsubj,size(ptime,2));

cfg                 = [];
cfg.avgoverchan     = 'yes';
cfg.frequency       = [8 12];
cfg.avgoverfreq     = 'yes';
sensorL = [121 122 139 140 123 124 127 128 129 130 141 142 143 144 147 148 149 150 155 156];
sensorR = [169 170 187 188 185 186 191 192 189 190 173 174 179 180 175 176 151 152 153 154];
sensorL = sensorL([2:2:numel(sensorL)])/2;
sensorR = sensorR([2:2:numel(sensorR)])/2;
cfg.channel         = sensorL;
for i=1:nsubj
    tfr_dataHiLo_left{i} = ft_selectdata(cfg, tfr_dataHiLo{i});
    tfr_dataHiLo_left{i}.label = {'MEG0102'}; %used as placeholder
    HiLo_left_sbj(i,:) = squeeze(tfr_dataHiLo_left{i}.powspctrm);
    tfr_dataLoHi_left{i} = ft_selectdata(cfg, tfr_dataLoHi{i});
    tfr_dataLoHi_left{i}.label = {'MEG0102'}; %used as placeholder
    LoHi_left_sbj(i,:) = squeeze(tfr_dataLoHi_left{i}.powspctrm);
end
cfg.channel         = sensorR;
for i=1:nsubj
    tfr_dataHiLo_right{i} = ft_selectdata(cfg, tfr_dataHiLo{i});
    tfr_dataHiLo_right{i}.label = {'MEG0102'}; %used as placeholder
    HiLo_right_sbj(i,:) = squeeze(tfr_dataHiLo_right{i}.powspctrm);
    tfr_dataLoHi_right{i} = ft_selectdata(cfg, tfr_dataLoHi{i});
    tfr_dataLoHi_right{i}.label = {'MEG0102'}; %used as placeholder
    LoHi_right_sbj(i,:) = squeeze(tfr_dataLoHi_right{i}.powspctrm);
end

% average across hemispheres
tfr_dataHi = cell(nsubj,1);
tfr_dataLo = cell(nsubj,1);
alphatime_Hi_sbj = nan(nsubj,size(ptime,2));
alphatime_Lo_sbj = nan(nsubj,size(ptime,2));
cfg             = [];
cfg.parameter   = 'powspctrm';
cfg.operation   = '(x1+x2)/2';
for i=1:nsubj
    tfr_dataHi{i} = ft_math(cfg, tfr_dataHiLo_right{i}, tfr_dataLoHi_left{i});
    alphatime_Hi_sbj(i,:) = squeeze(tfr_dataHi{i}.powspctrm);
    tfr_dataLo{i} = ft_math(cfg, tfr_dataHiLo_left{i}, tfr_dataLoHi_right{i});
    alphatime_Lo_sbj(i,:) = squeeze(tfr_dataLo{i}.powspctrm);
end

% grand average
cfg         = [];
tfr_dataHi_ga = ft_freqgrandaverage(cfg, tfr_dataHi{:});
tfr_dataLo_ga = ft_freqgrandaverage(cfg, tfr_dataLo{:});
alphatime_Hi_ga = squeeze(tfr_dataHi_ga.powspctrm);
alphatime_Lo_ga = squeeze(tfr_dataLo_ga.powspctrm);

% standard error
tfr_dataHi_ga.sem = std(alphatime_Hi_sbj)/sqrt(nsubj);
tfr_dataLo_ga.sem = std(alphatime_Lo_sbj)/sqrt(nsubj);


% Plot for Fig. 5 of the paper
tfr_dataHi_ga.mask = stat.mask; % adding mask to data
tfr_dataLo_ga.mask = stat.mask; % adding mask to data

cfg               = [];
cfg.maskparameter = 'mask';
cfg.maskstyle     = 'box';
cfg.maskfacealpha = 0.5; % transparency of mask
if strcmp(baselinecor, 'yes')
    cfg.ylim        = [1  2.1];
else
    cfg.ylim        = [.15e-23 1.2e-23];
end
figure(); ft_singleplotER(cfg, tfr_dataHi_ga, tfr_dataLo_ga); hold on;

ptime = tfr_dataHi_ga.time;
x_HiLo = ptime(1:60);
y1 = squeeze(tfr_dataHi_ga.powspctrm - tfr_dataHi_ga.sem)'; y1 = y1(1:60);
y2 = squeeze(tfr_dataHi_ga.powspctrm + tfr_dataHi_ga.sem)'; y2 = y2(1:60);
fill([x_HiLo fliplr(x_HiLo)], [y1 fliplr(y2)], 'b', 'LineStyle','none'); %plot standard error
alpha(0.10)
y1 = squeeze(tfr_dataLo_ga.powspctrm - tfr_dataLo_ga.sem)'; y1 = y1(1:60);
y2 = squeeze(tfr_dataLo_ga.powspctrm + tfr_dataLo_ga.sem)'; y2 = y2(1:60);
fill([x_HiLo fliplr(x_HiLo)], [y1 fliplr(y2)], 'r', 'LineStyle','none'); hold off %plot standard error
alpha(0.10)
xlabel('time'); ylabel('Alpha power');
xline(0, '--'); xline(1.5, '--'); 
yline(0, '--'); title('');
legend('high','low'); legend boxoff


%% Alpha-powr analysis using z-score baseline correction (Wang et al. JoCN 2019) 

clear

% set params
folder = 'Q:/MEG_data/';
filelist = dir(fullfile(folder, '**', 'tfr_kt3.mat')); %eye-tracking controled

filelistHiLo = filelist([1 5 3 8  9 13 12 16 17 20 21 24 25 29 27 34 32 36 37 40]); %HIGH-LOW condition
filelistLoHi = filelist([7 2 6 4 11 10 15 14 19 18 23 22 28 26 31 30 35 33 39 38]); %LOW-HIGH condition
nsubj = size(filelistLoHi,1);

% load data
tfr_dataHiLo = cell(nsubj,1);
tfr_dataLoHi = cell(nsubj,1);
for i=1:nsubj
    fprintf(1,'#%d\n',i); %print subject number
    
    % high-low cond
    tfr_dataHiLo{i} = load([filelistHiLo(i).folder '\' filelistHiLo(i).name]);
    % select frequency range
    cfg             = [];
    cfg.frequency   = [8 12];
    tfr_dataHiLo{i} = ft_selectdata(cfg, tfr_dataHiLo{i}.tfr_all);
    % combine planar gradiometers
    cfg                 = [];
    cfg.method          = 'sum';
    tfr_dataHiLo{i} = ft_combineplanar(cfg,tfr_dataHiLo{i});

    % low-high cond
    tfr_dataLoHi{i} = load([filelistLoHi(i).folder '\' filelistLoHi(i).name]);
    % select frequency range
    cfg             = [];
    cfg.frequency   = [8 12];
    tfr_dataLoHi{i} = ft_selectdata(cfg, tfr_dataLoHi{i}.tfr_all);
    % combine planar gradiometers
    cfg                 = [];
    cfg.method          = 'sum';
    tfr_dataLoHi{i} = ft_combineplanar(cfg,tfr_dataLoHi{i});
end

% selected sensors matched by pairs across hemifields
left_sensors  = [121 122 139 140 123 124 127 128 129 130 141 142 143 144 145 146 147 148 155 156];
right_sensors = [169 170 187 188 185 186 191 192 189 190 173 174 179 180 177 178 175 176 153 154];
left_sensors = left_sensors([2:2:numel(left_sensors)])/2;  %combine sens
right_sensors = right_sensors([2:2:numel(right_sensors)])/2;  %combine sens
cont_high_HiLo = cell(nsubj,1);
ipsi_high_HiLo = cell(nsubj,1);
cont_high_LoHi = cell(nsubj,1);
ipsi_high_LoHi = cell(nsubj,1);
for i = 1:nsubj
    % select right sensors
    cfg = [];
    cfg.channel = right_sensors;
    cont_high_HiLo{i} = ft_selectdata(cfg, tfr_dataHiLo{i});
    ipsi_high_LoHi{i} = ft_selectdata(cfg, tfr_dataLoHi{i});
    % select left sensors
    cfg = [];
    cfg.channel = left_sensors;
    ipsi_high_HiLo{i} = ft_selectdata(cfg, tfr_dataHiLo{i});
    ipsi_high_HiLo{i}.label = cont_high_HiLo{i}.label;  % renamed sensor labels for consistency 
    cont_high_LoHi{i} = ft_selectdata(cfg, tfr_dataLoHi{i});
    cont_high_LoHi{i}.label = cont_high_HiLo{i}.label;  % renamed sensor labels for consistency 
end

% compute difference between contra and ipsi
cont_ipsi_dif_HiLo = cell(nsubj,1);
cont_ipsi_dif_LoHi = cell(nsubj,1);
for i = 1:nsubj
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'subtract';
    cont_ipsi_dif_HiLo{i} = ft_math(cfg, cont_high_HiLo{i}, ipsi_high_HiLo{i});
    cont_ipsi_dif_LoHi{i} = ft_math(cfg, cont_high_LoHi{i}, ipsi_high_LoHi{i});
end

% run permutation
av_HiLo = zeros(1000,nsubj,10,5,61);  % nperm * nsubj * nsens * nfreq * ntime
av_LoHi = zeros(1000,nsubj,10,5,61);  % nperm * nsubj * nsens * nfreq * ntime
for i = 1:nsubj
    fprintf('#%d\n',i)
    for p = 1:1000
        % get data
        cont_data_HiLo = cont_high_HiLo{i}.powspctrm;
        ipsi_data_HiLo = ipsi_high_HiLo{i}.powspctrm;
        cont_data_LoHi = cont_high_LoHi{i}.powspctrm;
        ipsi_data_LoHi = ipsi_high_LoHi{i}.powspctrm;
        % get nr of trials
        ntrials_HiLo = size(cont_data_HiLo,1);
        ntrials_LoHi = size(cont_data_LoHi,1);
        % shuffle conditions across trials
        idx_HiLo = randperm(ntrials_HiLo,fix(ntrials_HiLo/2));
        x_HiLo = cat(1,cont_data_HiLo(idx_HiLo,:,:,:), ipsi_data_HiLo(setdiff(1:end,idx_HiLo),:,:,:));
        y_HiLo = cat(1,ipsi_data_HiLo(idx_HiLo,:,:,:), cont_data_HiLo(setdiff(1:end,idx_HiLo),:,:,:));
        idx_LoHi = randperm(ntrials_LoHi,fix(ntrials_LoHi/2));
        x_LoHi = cat(1,cont_data_LoHi(idx_LoHi,:,:,:), ipsi_data_LoHi(setdiff(1:end,idx_LoHi),:,:,:));
        y_LoHi = cat(1,ipsi_data_LoHi(idx_LoHi,:,:,:), cont_data_LoHi(setdiff(1:end,idx_LoHi),:,:,:));
        % compute difference between surrogate contra and ipsi
        d_HiLo = x_HiLo-y_HiLo;
        d_LoHi = x_LoHi-y_LoHi;
        % compute averaged difference
        av_HiLo(p,i,:,:,:) = squeeze(mean(d_HiLo,1));
        av_LoHi(p,i,:,:,:) = squeeze(mean(d_LoHi,1));
    end
end

% transform raw contra-ipsi differences into z-scores
cont_ipsi_dif_HiLo_z = cont_ipsi_dif_HiLo;
cont_ipsi_dif_LoHi_z = cont_ipsi_dif_LoHi;
for i = 1:nsubj
    x_HiLo = squeeze(cont_ipsi_dif_HiLo{i}.powspctrm);
    m_HiLo = squeeze(mean(av_HiLo(:,i,:,:,:),1));
    s_HiLo = squeeze(std(av_HiLo(:,i,:,:,:),0,1));
    for t = size(x_HiLo,1)
        cont_ipsi_dif_HiLo_z{i}.powspctrm(t,:,:,:) = (squeeze(x_HiLo(t,:,:,:))-m_HiLo)./s_HiLo;
    end
    x_LoHi = squeeze(cont_ipsi_dif_LoHi{i}.powspctrm);
    m_LoHi = squeeze(mean(av_LoHi(:,i,:,:,:),1));
    s_LoHi = squeeze(std(av_LoHi(:,i,:,:,:),0,1));
    for t = size(x_LoHi,1)
        cont_ipsi_dif_LoHi_z{i}.powspctrm(t,:,:,:) = (squeeze(x_LoHi(t,:,:,:))-m_LoHi)./s_LoHi;
    end
end

% average across trials, sensors and frequencies
cont_ipsi_dif_HiLo_z_av = cell(nsubj,1);
cont_ipsi_dif_LoHi_z_av = cell(nsubj,1);
cfg = [];
cfg.avgoverrpt = 'yes';
cfg.avgoverchan = 'yes';
cfg.avgoverfreq = 'yes';
for i = 1:nsubj
    cont_ipsi_dif_HiLo_z_av{i} = ft_selectdata(cfg, cont_ipsi_dif_HiLo_z{i});
    cont_ipsi_dif_LoHi_z_av{i} = ft_selectdata(cfg, cont_ipsi_dif_LoHi_z{i});
end

% grand average within sessions
cfg         = [];
cont_ipsi_dif_HiLo_ga = ft_freqgrandaverage(cfg, cont_ipsi_dif_HiLo_z_av{:});
cont_ipsi_dif_LoHi_ga = ft_freqgrandaverage(cfg, cont_ipsi_dif_LoHi_z_av{:});

% average across sessions
cont_ipsi_dif_z_av = cont_ipsi_dif_HiLo_z_av;
for i = 1:nsubj
    cont_ipsi_dif_z_av{i}.powspctrm = mean([ ...
        cont_ipsi_dif_HiLo_z_av{i}.powspctrm, ...
        cont_ipsi_dif_LoHi_z_av{i}.powspctrm],2);
end

% grand average across sessions
cont_ipsi_dif_ga = ft_freqgrandaverage([], cont_ipsi_dif_z_av{:});

% average within the placeholder time window and get individual data
cont_ipsi_dif_z_av_sbj = nan(nsubj,1);
for i = 1:nsubj
    cfg             = [];
    cfg.latency = [0 1.5];
    cfg.avgovertime = 'yes';
    temp = ft_selectdata(cfg, cont_ipsi_dif_z_av{i});
    cont_ipsi_dif_z_av_sbj(i) = temp.powspctrm;
end


% plot (not in the paper)
figure
violinplot(cont_ipsi_dif_z_av_sbj, [], 'ShowMean', true); hold on
yline(0, '--'); hold off
xlim([0.5,1.5]); title('high probability: contra minus ipsi');
ylim([-.02, .02]); ylabel('Alpha power (z-score)');
