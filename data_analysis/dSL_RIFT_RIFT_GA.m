%% RIFT cross-correlation analysis

clear
ft_defaults

% load xcorr data
folder = 'Q:/MEG_data/';
filelistL = dir(fullfile(folder, '**', 'xcorL3.mat'));
filelistR = dir(fullfile(folder, '**', 'xcorR3.mat'));

filelistHiL = filelistL([1 5 8  9 13 12 16 17 20 21 24 25 29 27 34 32 36 37 40]); %HIGH-LOW condition (sbj3 excluded)
filelistLoR = filelistR([1 5 8  9 13 12 16 17 20 21 24 25 29 27 34 32 36 37 40]); %HIGH-LOW condition (sbj3 excluded)
filelistLoL = filelistL([7 2 4 11 10 15 14 19 18 23 22 28 26 31 30 35 33 39 38]); %LOW-HIGH condition (sbj3 excluded)
filelistHiR = filelistR([7 2 4 11 10 15 14 19 18 23 22 28 26 31 30 35 33 39 38]); %LOW-HIGH condition (sbj3 excluded)

xcorHiL_data = zeros(204,401,size(filelistHiL,1));
xcorLoR_data = zeros(204,401,size(filelistLoR,1));
xcorLoL_data = zeros(204,401,size(filelistLoL,1));
xcorHiR_data = zeros(204,401,size(filelistHiR,1));
for i = 1:size(filelistHiL,1)
    xcorHiL_temp = load([filelistHiL(i).folder '\' filelistHiL(i).name]);
    xcorHiL_data(:,:,i) = cell2mat(struct2cell(xcorHiL_temp));
    xcorLoR_temp = load([filelistLoR(i).folder '\' filelistLoR(i).name]);
    xcorLoR_data(:,:,i) = cell2mat(struct2cell(xcorLoR_temp));
    xcorLoL_temp = load([filelistLoL(i).folder '\' filelistLoL(i).name]);
    xcorLoL_data(:,:,i) = cell2mat(struct2cell(xcorLoL_temp));
    xcorHiR_temp = load([filelistHiR(i).folder '\' filelistHiR(i).name]);
    xcorHiR_data(:,:,i) = cell2mat(struct2cell(xcorHiR_temp));
end

% compute hilbert transform
xcorHiL_hilb = zeros(401,204,size(filelistHiL,1));
xcorLoR_hilb = zeros(401,204,size(filelistHiL,1));
xcorLoL_hilb = zeros(401,204,size(filelistHiL,1));
xcorHiR_hilb = zeros(401,204,size(filelistHiL,1));
for i = 1:size(filelistHiL,1)
    xcorHiL_hilb(:,:,i) = abs(hilbert(xcorHiL_data(:,:,i)'));
    xcorLoR_hilb(:,:,i) = abs(hilbert(xcorLoR_data(:,:,i)'));
    xcorLoL_hilb(:,:,i) = abs(hilbert(xcorLoL_data(:,:,i)'));
    xcorHiR_hilb(:,:,i) = abs(hilbert(xcorHiR_data(:,:,i)'));
end

% normalize
xcorHiL_hilb_n = zeros(401,204,size(filelistHiL,1));
xcorLoR_hilb_n = zeros(401,204,size(filelistHiL,1));
xcorLoL_hilb_n = zeros(401,204,size(filelistHiL,1));
xcorHiR_hilb_n = zeros(401,204,size(filelistHiL,1));
for i = 1:size(filelistHiL,1)
    minHi = min([xcorHiL_hilb(:,:,i), xcorLoR_hilb(:,:,i)],[],'all'); maxHL = max([xcorHiL_hilb(:,:,i), xcorLoR_hilb(:,:,i)],[],'all');
    xcorHiL_hilb_n(:,:,i) = (xcorHiL_hilb(:,:,i) - minHi) / (maxHL - minHi);
    xcorLoR_hilb_n(:,:,i) = (xcorLoR_hilb(:,:,i) - minHi) / (maxHL - minHi);
    minHi = min([xcorLoL_hilb(:,:,i), xcorHiR_hilb(:,:,i)],[],'all'); maxHL = max([xcorLoL_hilb(:,:,i), xcorHiR_hilb(:,:,i)],[],'all');
    xcorLoL_hilb_n(:,:,i) = (xcorLoL_hilb(:,:,i) - minHi) / (maxHL - minHi);
    xcorHiR_hilb_n(:,:,i) = (xcorHiR_hilb(:,:,i) - minHi) / (maxHL - minHi);
end

% get peak value
xcorHiL_peak = max(xcorHiL_hilb_n,[],1);
xcorLoR_peak = max(xcorLoR_hilb_n,[],1);
xcorLoL_peak = max(xcorLoL_hilb_n,[],1);
xcorHiR_peak = max(xcorHiR_hilb_n,[],1);

% compute difference between stimuli
xcorHLdiff = zeros(1,204,size(filelistHiL,1));
xcorLHdiff = zeros(1,204,size(filelistHiL,1));
for i = 1:size(filelistHiL,1)
    xcorHLdiff(:,:,i) = (xcorHiL_peak(:,:,i) - xcorLoR_peak(:,:,i));% ./ (xcorHiL_peak(:,:,i) + xcorLoR_peak(:,:,i));
    xcorLHdiff(:,:,i) = (xcorLoL_peak(:,:,i) - xcorHiR_peak(:,:,i));% ./ (xcorLoL_peak(:,:,i) + xcorHiR_peak(:,:,i));
end

% select the best sensor for each stimulation side (left/negative and right/positive)
sens_nr = 1;
xcorHLdiff_Rtop = zeros(1,sens_nr,size(filelistHiL,1));
xcorHLdiff_Ltop = zeros(1,sens_nr,size(filelistHiL,1));
xcorHLdiff_sel = zeros(1,204,size(filelistHiL,1));
xcorLHdiff_Rtop = zeros(1,sens_nr,size(filelistHiL,1));
xcorLHdiff_Ltop = zeros(1,sens_nr,size(filelistHiL,1));
xcorLHdiff_sel = zeros(1,204,size(filelistHiL,1));
for i = 1:size(filelistHiL,1)
    xcorHLdiff_Rtop(:,:,i) = mink(xcorHLdiff(:,:,i),sens_nr); xcorHLdiff_Ltop(:,:,i) = maxk(xcorHLdiff(:,:,i),sens_nr);
    xcorHLdiff_sel(:,:,i) = xcorHLdiff(:,:,i);
    xcorHLdiff_sel(:,~ismember(xcorHLdiff_sel(:,:,i), [xcorHLdiff_Ltop(:,:,i), xcorHLdiff_Rtop(:,:,i)]),i) = 0;
    xcorLHdiff_Rtop(:,:,i) = mink(xcorLHdiff(:,:,i),sens_nr); xcorLHdiff_Ltop(:,:,i) = maxk(xcorLHdiff(:,:,i),sens_nr);
    xcorLHdiff_sel(:,:,i) = xcorLHdiff(:,:,i);
    xcorLHdiff_sel(:,~ismember(xcorLHdiff_sel(:,:,i), [xcorLHdiff_Ltop(:,:,i), xcorLHdiff_Rtop(:,:,i)]),i) = 0;
end

% compute mean
xcorHLdiff_Rtop_sbj = mean(xcorHLdiff_Rtop,2);
xcorHLdiff_Rtop_mean = mean(xcorHLdiff_Rtop_sbj,3);
xcorHLdiff_Rtop_sem = std(xcorHLdiff_Rtop_sbj,[],3)/sqrt(size(filelistHiL,1));
xcorHLdiff_Ltop_sbj = mean(xcorHLdiff_Ltop,2);
xcorHLdiff_Ltop_mean = mean(xcorHLdiff_Ltop_sbj,3);
xcorHLdiff_Ltop_sem = std(xcorHLdiff_Ltop_sbj,[],3)/sqrt(size(filelistHiL,1));

xcorLHdiff_Rtop_sbj = mean(xcorLHdiff_Rtop,2);
xcorLHdiff_Rtop_mean = mean(xcorLHdiff_Rtop_sbj,3);
xcorLHdiff_Rtop_sem = std(xcorLHdiff_Rtop_sbj,[],3)/sqrt(size(filelistHiL,1));
xcorLHdiff_Ltop_sbj = mean(xcorLHdiff_Ltop,2);
xcorLHdiff_Ltop_mean = mean(xcorLHdiff_Ltop_sbj,3);
xcorLHdiff_Ltop_sem = std(xcorLHdiff_Ltop_sbj,[],3)/sqrt(size(filelistHiL,1));


% Plot for Fig. 4 of the paper [part 1]
addpath(['...\Violinplot'])
% plot xcor by condition (hemi x probability)
f1 = figure('Units', 'centimeters', 'Position', [5,5,8,8]);
subplot(221)
xcorHLdiff_sbj = [squeeze(xcorHLdiff_Ltop_sbj), squeeze(xcorLHdiff_Ltop_sbj),...
    squeeze(-xcorLHdiff_Rtop_sbj), squeeze(-xcorHLdiff_Rtop_sbj)];
vp = violinplot(xcorHLdiff_sbj, [], 'ShowMean', true);
xlim([0.5,4.5]); xlabel('Distractor probability'); xticklabels({'high', 'low', 'high', 'low'});
ylim([0 1.1]); ylabel('Frequency tagging response');
legend('left','right','FontSize',6); legend boxoff
ax = gca;
ax.FontSize = 6;
% set colours
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

% compute difference between conditions
xcorH_sbj = mean([xcorHLdiff_Ltop_sbj, -xcorLHdiff_Rtop_sbj],2);
xcorH_sbj = squeeze(xcorH_sbj);
xcorL_sbj = mean([-xcorHLdiff_Rtop_sbj, xcorLHdiff_Ltop_sbj],2);
xcorL_sbj = squeeze(xcorL_sbj);

%  compute statistics
addpath('C:\Users\ferranto\DataShare\UoB\Scripts\Bayes_factor')
addpath('C:\Users\ferranto\DataShare\UoB\Scripts\Cohens_d')
[~,p,~,stats] = ttest(xcorH_sbj, xcorL_sbj) %t-test
d_HiLo = computeCohen_d(xcorH_sbj, xcorL_sbj, 'paired')
% [bf10,~] = bf.ttest(xcorH_sbj, xcorL_sbj)


% Plot for Fig. 4 of the paper [part 2]
% plot xcor probability effect (across hemi)
figure(f1)
subplot(222)
xcor_sbj = squeeze(xcorL_sbj) - squeeze(xcorH_sbj);
violinplot(xcor_sbj, [], 'ShowMean', true); hold on
yline(0, '--'); hold off
xlim([0.5,1.5]); xticklabels('low-high');
ylim([-.15 .35]); ylabel('SL effect (norm corr)');
ax = gca;
ax.FontSize = 6;



%% RIFT cross-correlation analysis splitting the session data into 3 parts

clear

% set folder and select files
folder = 'Q:/MEG_data/';
filelistL = dir(fullfile(folder, '**', 'xcorLkt3.mat')); %eye-tracking controled
filelistR = dir(fullfile(folder, '**', 'xcorRkt3.mat')); %eye-tracking controled

filelistHiL = filelistL([1 5 8  9 13 12 16 17 20 21 24 25 29 27 34 32 36 37 40]); %HIGH-LOW condition (sbj3 excluded)
filelistLoR = filelistR([1 5 8  9 13 12 16 17 20 21 24 25 29 27 34 32 36 37 40]); %HIGH-LOW condition (sbj3 excluded)
filelistLoL = filelistL([7 2 4 11 10 15 14 19 18 23 22 28 26 31 30 35 33 39 38]); %LOW-HIGH condition (sbj3 excluded)
filelistHiR = filelistR([7 2 4 11 10 15 14 19 18 23 22 28 26 31 30 35 33 39 38]); %LOW-HIGH condition (sbj3 excluded)
nsubj = size(filelistHiL,1);

% load data
xcorHiL_data = cell(nsubj,1);
xcorLoR_data = cell(nsubj,1);
xcorLoL_data = cell(nsubj,1);
xcorHiR_data = cell(nsubj,1);
for i = 1:nsubj
    fprintf(1,'#%d\n',i);
    xcorHiL_data{i} = load([filelistHiL(i).folder '\' filelistHiL(i).name]);
    xcorLoR_data{i} = load([filelistLoR(i).folder '\' filelistLoR(i).name]);
    xcorLoL_data{i} = load([filelistLoL(i).folder '\' filelistLoL(i).name]);
    xcorHiR_data{i} = load([filelistHiR(i).folder '\' filelistHiR(i).name]);
end

% split trials into 3 epochs (1=early, 2=middle and 3=late) and average
xcorHiL_data_split = cell(nsubj,3);
xcorLoR_data_split = cell(nsubj,3);
xcorLoL_data_split = cell(nsubj,3);
xcorHiR_data_split = cell(nsubj,3);
for i = 1:nsubj
    ntrialsHiL_LoR = size(xcorHiL_data{i}.xcorL,3);
    qsHiL_LoR = quantile(1:ntrialsHiL_LoR,[0, 1/3, 2/3, 1]);
    ntrialsLoL_HiR = size(xcorLoL_data{i}.xcorL,3);
    qsLoL_HiR = quantile(1:ntrialsLoL_HiR,[0, 1/3, 2/3, 1]);
    for j = 1:3
        t_idxHiL_LoR = [round(qsHiL_LoR(j)):round(qsHiL_LoR(j+1))-1];
        xcorHiL_data_split{i,j} = mean(xcorHiL_data{i}.xcorL(:,:,t_idxHiL_LoR),3);
        xcorLoR_data_split{i,j} = mean(xcorLoR_data{i}.xcorR(:,:,t_idxHiL_LoR),3);
        t_idxLoL_HiR = [round(qsLoL_HiR(j)):round(qsLoL_HiR(j+1))-1];
        xcorLoL_data_split{i,j} = mean(xcorLoL_data{i}.xcorL(:,:,t_idxLoL_HiR),3);
        xcorHiR_data_split{i,j} = mean(xcorHiR_data{i}.xcorR(:,:,t_idxLoL_HiR),3);
    end
end

% compute hilbert transform
xcorHiL_hilb = cell(nsubj,3);
xcorLoR_hilb = cell(nsubj,3);
xcorLoL_hilb = cell(nsubj,3);
xcorHiR_hilb = cell(nsubj,3);
for i = 1:nsubj
    for j = 1:3
        xcorHiL_hilb{i,j} = abs(hilbert(xcorHiL_data_split{i,j}'));
        xcorLoR_hilb{i,j} = abs(hilbert(xcorLoR_data_split{i,j}'));
        xcorLoL_hilb{i,j} = abs(hilbert(xcorLoL_data_split{i,j}'));
        xcorHiR_hilb{i,j} = abs(hilbert(xcorHiR_data_split{i,j}'));
    end
end

% get peak value
xcorHiL_peak = cell(nsubj,3);
xcorLoR_peak = cell(nsubj,3);
xcorLoL_peak = cell(nsubj,3);
xcorHiR_peak = cell(nsubj,3);
for i = 1:nsubj
    for j = 1:3
        xcorHiL_peak{i,j} = max(xcorHiL_hilb{i,j},[],1);
        xcorLoR_peak{i,j} = max(xcorLoR_hilb{i,j},[],1);
        xcorLoL_peak{i,j} = max(xcorLoL_hilb{i,j},[],1);
        xcorHiR_peak{i,j} = max(xcorHiR_hilb{i,j},[],1);
    end
end

% compute difference between stimuli
xcorHLdiff = cell(nsubj,3);
xcorLHdiff = cell(nsubj,3);
for i = 1:nsubj
    for j = 1:3
        xcorHLdiff{i,j} = (xcorHiL_peak{i,j} - xcorLoR_peak{i,j});  % ./ (xcorHiL_peak(:,:,i) + xcorLoR_peak(:,:,i));
        xcorLHdiff{i,j} = (xcorLoL_peak{i,j} - xcorHiR_peak{i,j});  % ./ (xcorLoL_peak(:,:,i) + xcorHiR_peak(:,:,i));
    end
end

% select top sensors of each stimulus side (left/negative and right/positive)
sens_nr = 3;
xcorHLdiff_Rtop = cell(nsubj,3);
xcorHLdiff_Ltop = cell(nsubj,3);
xcorLHdiff_Rtop = cell(nsubj,3);
xcorLHdiff_Ltop = cell(nsubj,3);
for i = 1:nsubj
    for j = 1:3
        xcorHLdiff_Rtop{i,j} = mink(xcorHLdiff{i,j},sens_nr); 
        xcorHLdiff_Ltop{i,j} = maxk(xcorHLdiff{i,j},sens_nr);
        xcorLHdiff_Rtop{i,j} = mink(xcorLHdiff{i,j},sens_nr); 
        xcorLHdiff_Ltop{i,j} = maxk(xcorLHdiff{i,j},sens_nr);
    end
end

% average over sensors
xcorHLdiff_Rtop_sbj = nan(nsubj,3);
xcorHLdiff_Ltop_sbj = nan(nsubj,3);
xcorLHdiff_Rtop_sbj = nan(nsubj,3);
xcorLHdiff_Ltop_sbj = nan(nsubj,3);
for j = 1:3
    xcorHLdiff_Rtop_sbj(:,j) = mean(cell2mat(xcorHLdiff_Rtop(:,j)),2);
    xcorHLdiff_Ltop_sbj(:,j) = mean(cell2mat(xcorHLdiff_Ltop(:,j)),2);
    xcorLHdiff_Rtop_sbj(:,j) = mean(cell2mat(xcorLHdiff_Rtop(:,j)),2);
    xcorLHdiff_Ltop_sbj(:,j) = mean(cell2mat(xcorLHdiff_Ltop(:,j)),2);
end

% compute difference
xcorH_sbj = mean(cat(3,abs(xcorHLdiff_Ltop_sbj), abs(xcorLHdiff_Rtop_sbj)),3);
xcorL_sbj = mean(cat(3,abs(xcorHLdiff_Rtop_sbj), abs(xcorLHdiff_Ltop_sbj)),3);


% plot (not in the paper)
addpath(['...\Scripts\Violinplot'])
% plot xcor by condition (hemi x probability)
f1 = figure('Units', 'centimeters', 'Position', [5,5,8,8]);
for j = 1:3
    subplot(2,3,j)
    xcorHLdiff_sbj = [squeeze(xcorHLdiff_Ltop_sbj(:,j)), squeeze(xcorLHdiff_Ltop_sbj(:,j)),...
        squeeze(-xcorLHdiff_Rtop_sbj(:,j)), squeeze(-xcorHLdiff_Rtop_sbj(:,j))];
    vp = violinplot(xcorHLdiff_sbj, [], 'ShowMean', true);
    xlim([0.5,4.5]); xlabel('Distractor probability'); xticklabels({'high', 'low', 'high', 'low'});
    ylabel('Frequency tagging response');
    legend('left','right','FontSize',6); legend boxoff
    ax = gca;
    ax.FontSize = 6;
    % set colours
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
    legend(hleglines, {'Session 1', 'Session 2'},'Location','northeast')
    legend('boxoff')
end
% plot xcor probability effect (across hemi)
figure(f1)
for j = 1:3
    subplot(2,3,j+3)
    xcor_sbj = squeeze(xcorL_sbj(:,j)) - squeeze(xcorH_sbj(:,j));
    violinplot(xcor_sbj, [], 'ShowMean', true); hold on
    ylim([-.02,.03]); yline(0, '--'); hold off
    xlim([0.5,1.5]); xticklabels('low-high');
    ylabel('SL effect (norm corr)');
    ax = gca;
    ax.FontSize = 6;
end


%% Correlation behaviour-bft-alpha --- load data first ----

beh_data = sbj_av_drt(:,[1:2,4:20]); %RTs
beh_data_diff = beh_data(2,:)' - beh_data(1,:)';

bft_data_diff = squeeze(xcorL_sbj - xcorH_sbj);

% compute correlation behaviour-rift
[rho,pval] = corr(beh_data_diff,bft_data_diff,'type','pearson')


% plot (not in the paper)
figure(); plot(beh_data_diff,bft_data_diff,'xk','MarkerSize',12); hold on; hl = lsline; set(hl,'linewidth',2);
xlabel('SL effect on RTs'); ylabel('SL effect on RIFT')
