
clearvars
clc

%% Load data

cd('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data'); %ferranto
filename = dir('*.txt');
filename(contains({filename.name}, '1_B64D_s2p1.txt')) = []; %remove '1_B64D_s2p1.txt'
startRow = 2;
endRow = 361;
sbj_nr = 20;

for i=1:numel(filename)
    temp = importfile(filename(i).name, startRow, endRow);
    if i == 1
        data = temp;
    else
        data = [data;temp];
    end
end


%% New variables

% Create distractor presence variable
data.distractor_pres = data.distractor_loc > 0;

% Create distractor-location probability variable
data.distractor_prob = repmat('a',size(data,1),1);
for i=1:size(data,1)
    if mod(data.subject_nr(i),2) == 1
        if data.session_count(i) < 2
            if ismember(data.distractor_loc(i),[1,2])
                data.distractor_prob(i) = 'h';
            elseif ismember(data.distractor_loc(i),[3,4])
                data.distractor_prob(i) = 'l';
            end
        else
            if ismember(data.distractor_loc(i),[1,2])
                data.distractor_prob(i) = 'l';
            elseif ismember(data.distractor_loc(i),[3,4])
                data.distractor_prob(i) = 'h';
            end            
        end
    else
        if data.session_count(i) < 2
            if ismember(data.distractor_loc(i),[1,2])
                data.distractor_prob(i) = 'l';
            elseif ismember(data.distractor_loc(i),[3,4])
                data.distractor_prob(i) = 'h';
            end
        else
            if ismember(data.distractor_loc(i),[1,2])
                data.distractor_prob(i) = 'h';
            elseif ismember(data.distractor_loc(i),[3,4])
                data.distractor_prob(i) = 'l';
            end
        end
    end
end

% Create target-location probability variable
data.target_prob = repmat('x',size(data,1),1);
for i=1:size(data,1)
    if mod(data.subject_nr(i),2) == 1
        if data.session_count(i) < 2
            if ismember(data.target_loc(i),[1,2])
                data.target_prob(i) = 'h';
            elseif ismember(data.target_loc(i),[3,4])
                data.target_prob(i) = 'l';
            end
        else
            if ismember(data.target_loc(i),[1,2])
                data.target_prob(i) = 'l';
            elseif ismember(data.target_loc(i),[3,4])
                data.target_prob(i) = 'h';
            end            
        end
    else
        if data.session_count(i) < 2
            if ismember(data.target_loc(i),[1,2])
                data.target_prob(i) = 'l';
            elseif ismember(data.target_loc(i),[3,4])
                data.target_prob(i) = 'h';
            end
        else
            if ismember(data.target_loc(i),[1,2])
                data.target_prob(i) = 'h';
            elseif ismember(data.target_loc(i),[3,4])
                data.target_prob(i) = 'l';
            end
        end
    end
end

% Create frequency-tagged distractor location variable
data.tagged_disloc = repmat('a',size(data,1),1);
for i=1:size(data,1)
    if ismember(data.distractor_loc(i),[2,4])
        data.tagged_disloc(i) = 'y';
    elseif ismember(data.distractor_loc(i),[1,3])
        data.tagged_disloc(i) = 'n';
    end
end

% Create frequency-tagged target location variable
data.tagged_tarloc = repmat('x',size(data,1),1);
for i=1:size(data,1)
    if ismember(data.target_loc(i),[2,4])
        data.tagged_tarloc(i) = 'y';
    elseif ismember(data.target_loc(i),[1,3])
        data.tagged_tarloc(i) = 'n';
    end
end

% Create inter-trial distractor-location repetition variable
data.repeat_disloc = repmat('n',size(data,1),1);
for n=unique(data.subject_nr)' %loop through subjects
    for s=1:2 %loop through sessions
        data_temp = data(data.subject_nr == n & data.session_count == s,:);
        for i=2:size(data_temp,1)
            if data_temp.distractor_loc(i) == 0 || data_temp.distractor_loc(i) ~= data_temp.distractor_loc(i-1)
                data_temp.repeat_disloc(i) = 'n';
            else
                data_temp.repeat_disloc(i) = 'y';
            end
        end
        data(data.subject_nr == n & data.session_count == s,'repeat_disloc') = table(data_temp.repeat_disloc);
    end
end


%% Tabulate

% [tbl,~,~,labels] = crosstab(data.subject_nr, data.distractor_loc, data.session_count) %distractor_loc
% [tbl,~,~,labels] = crosstab(data.subject_nr, data.distractor_prob, data.session_count) %distractor_prob
% [tbl,~,~,labels] = crosstab(data.distractor_loc, data.distractor_prob, data.session_count, data.subject_nr) %distractor_loc x distractor_prob

% [tbl,~,~,labels] = crosstab(data.subject_nr, data.target_loc, data.session_count) %target_loc
% [tbl,~,~,labels] = crosstab(data.subject_nr, data.target_loc, data.distractor_loc, data.session_count) %target_loc x distractor_loc


%% Remove initial 100 trials from all and uncorrect and slow response trials from RTs

% data_all = data; %data = data_all;
correct_trials = find(data.correct == 1);
nottoofast_trials = find(data.RT > 201);
nottooslow_trials = find(data.RT < mean(data.RT)+(2*std(data.RT)));
first100trials = find(~(data.trial_count <= 100 & data.part_count == 1));
good_trials = intersect(...
    first100trials, intersect(...
    correct_trials, intersect(...
    nottoofast_trials, nottooslow_trials)));

datart = data(good_trials,:);
data = data(first100trials,:);


%% Distractor location - RT analysis

% Mean reaction time
[means,stds,grp] = grpstats(datart.RT,...
    {datart.subject_nr, datart.session_count, datart.distractor_prob},...
    {'mean','std','gname'});
means2_rt = reshape(means,3,2,sbj_nr);
stds2 = reshape(stds,3,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),3,2,sbj_nr);

sbj_s1_rt = squeeze(means2_rt(:,1,:));
sbj_s2_rt = squeeze(means2_rt(:,2,:));
sbj_av_rt = squeeze(mean(means2_rt,2));

i_a = find(cell2mat(grp(1:3,3)) == 'a'); %find indexes of distractor probs
i_h = find(cell2mat(grp(1:3,3)) == 'h');
i_l = find(cell2mat(grp(1:3,3)) == 'l');

sbj_s1_drt = sbj_s1_rt([i_h,i_l],:) - sbj_s1_rt(i_a,:);
sbj_s2_drt = sbj_s2_rt([i_h,i_l],:) - sbj_s2_rt(i_a,:);
sbj_av_drt = sbj_av_rt([i_h,i_l],:) - sbj_av_rt(i_a,:);

i_h_d = 1;
i_l_d = 2;

group_s1_drt = mean(sbj_s1_drt,2);
group_s2_drt = mean(sbj_s2_drt,2);
group_av_drt = mean(sbj_av_drt,2);


% Plots for Fig. 3 of the paper [part 1]
addpath('...\Violinplot')
f = figure('Units', 'centimeters', 'Position', [5,5,8,8]);
% PLOT 1
% make first plot
subplot(221)
vp = violinplot([sbj_s1_drt', sbj_s2_drt'], [], 'ShowMean', true); hold on;
% set lalbels
xlim([0.5,4.5]); xlabel('Distarctor probability'); xticks(1:4); xticklabels({'high','low','high','low'});
ylim([10, 135]); ylabel('Distractor cost (ms)');
title('RTs');
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
% PLOT 2
% delta RTs between conditions (low-high)
sbj_av_drt_d = sbj_av_drt(i_l_d,:) - sbj_av_drt(i_h_d,:);
% make second plot
subplot(223)
violinplot(sbj_av_drt_d', [], 'ShowMean', true); hold on
yline(0, '--'); hold off
xlim([0.5,1.50]); xticklabels({'Low-High'});
ylim([-15, 50]); ylabel('SL effect (ms)');
ax = gca;
ax.FontSize = 6;
hold off


%% Distractor location - Accuracy analysis

% Mean accuracy
[means,stds,grp] = grpstats(data.correct,...
    {data.subject_nr, data.session_count, data.distractor_prob},...
    {'mean','std','gname'});
means2_ac = reshape(means,3,2,sbj_nr);
stds2 = reshape(stds,3,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),3,2,sbj_nr);

sbj_s1_acc = squeeze(means2_ac(:,1,:));
sbj_s2_acc = squeeze(means2_ac(:,2,:));
sbj_av_acc = squeeze(mean(means2_ac,2));

i_a = find(cell2mat(grp(1:3,3)) == 'a'); %find indexes of distractor probs
i_h = find(cell2mat(grp(1:3,3)) == 'h');
i_l = find(cell2mat(grp(1:3,3)) == 'l');

sbj_s1_dacc = sbj_s1_acc(i_a,:) - sbj_s1_acc([i_h,i_l],:);
sbj_s2_dacc = sbj_s2_acc(i_a,:) - sbj_s2_acc([i_h,i_l],:);
sbj_av_dacc = sbj_av_acc(i_a,:) - sbj_av_acc([i_h,i_l],:);

i_h_d = 1;
i_l_d = 2;

group_s1_dacc = mean(sbj_s1_dacc,2);
group_s2_dacc = mean(sbj_s2_dacc,2);
group_av_dacc = mean(sbj_av_dacc,2);


% Plot for Fig. 3 of the paper [part 2] 
figure(f)
% PLOT 1
% make first plot
subplot(222)
vp = violinplot([(sbj_s1_dacc*100)', (sbj_s2_dacc*100)'], [], 'ShowMean', true); hold on;
% set lalbels
xlim([0.5,4.5]); xlabel('Distarctor probability'); xticks(1:4); xticklabels({'high','low','high','low'});
ylim([-5, 20]); ylabel('Distractor cost (%acc)');
title('Accuracy');
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
% PLOT 2
% delta accuracy between conditions (high-low)
sbj_av_dacc_d = sbj_av_dacc(2,:) - sbj_av_dacc(1,:);
% make first plot
subplot(224)
violinplot((sbj_av_dacc_d*100)', [], 'ShowMean', true); hold on
yline(0, '--'); hold off
xlim([0.5,1.50]); xticklabels({'Low-High'});
ylim([-3, 10]); ylabel('SL effect (%acc)');
ax = gca;
ax.FontSize = 6;
hold off


%% Distractor location x Block - RT analysis

% exclude subject 1
datart2 = datart(datart.subject_nr ~= 1,:);
sbj_nr_minus1 = sbj_nr - 1;

% block_count values are reset at the beginning of part 2 within each
% experiment. This fix make them continues between the two part of the exp
datart2.block_count2 = zeros(size(datart2,1),1);
datart2.block_count2(datart2.part_count == 1) = datart2.block_count(datart2.part_count == 1);
datart2.block_count2(datart2.part_count == 2) = datart2.block_count(datart2.part_count == 2) + 30;

% each miniblock has 12 trials. to avoid having missing cells, here we
% average 3 miniblocks together
datart2.block_count3 = zeros(size(datart2,1),1);
datart2.block_count3 = round(datart2.block_count2 / 2);

% Mean reaction time
[means,stds,grp] = grpstats(datart2.RT,...
    {datart2.subject_nr, datart2.session_count, datart2.block_count3, datart2.distractor_prob},...
    {'mean','std','gname'});
means2 = reshape(means,3,30,2,sbj_nr_minus1);
stds2 = reshape(stds,3,30,2,sbj_nr_minus1);
grp2 = reshape(strcat(grp(:,1),'_',grp(:,2),'_',grp(:,3),'_',grp(:,4)),3,30,2,sbj_nr_minus1);

i_a = find(cell2mat(grp(1:3,4)) == 'a'); %find indexes of distractor probs
i_h = find(cell2mat(grp(1:3,4)) == 'h');
i_l = find(cell2mat(grp(1:3,4)) == 'l');

sbj_s1_rt = squeeze(means2(:,:,1,:));
sbj_s2_rt = squeeze(means2(:,:,2,:));
sbj_av_rt = squeeze(mean(means2,3));

group_s1_rt = mean(sbj_s1_rt,3);
group_s1_sem = std(sbj_s1_rt,[],3) / sqrt(sbj_nr_minus1);
group_s2_rt = mean(sbj_s2_rt,3);
group_s2_sem = std(sbj_s2_rt,[],3) / sqrt(sbj_nr_minus1);
group_av_rt = mean(sbj_av_rt,3);
group_av_sem = std(sbj_av_rt,[],3) / sqrt(sbj_nr_minus1);


% plot (not in the paper)
figure;
subplot(211);
errorbar(1:30,squeeze(group_s1_rt(i_a,:,:)), squeeze(group_s1_sem(i_a,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s1_rt(i_h,:,:)), squeeze(group_s1_sem(i_h,:,:)),'b');
errorbar(1:30,squeeze(group_s1_rt(i_l,:,:)), squeeze(group_s1_sem(i_l,:,:)),'g');
xlim([0.5, 30.5]); ylim([400, 670]); ylabel('RTs'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Session 1');
subplot(212);
errorbar(1:30,squeeze(group_s2_rt(i_a,:,:)), squeeze(group_s2_sem(i_a,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s2_rt(i_h,:,:)), squeeze(group_s2_sem(i_h,:,:)),'b');
errorbar(1:30,squeeze(group_s2_rt(i_l,:,:)), squeeze(group_s2_sem(i_l,:,:)),'g');
xlim([0.5, 30.5]); ylim([400, 670]); ylabel('RTs'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Session 2');
% subplot(313);
% errorbar(1:30,squeeze(group_av_rt(i_a,:,:)), squeeze(group_av_sem(i_a,:,:)),'r'); hold on;
% errorbar(1:30,squeeze(group_av_rt(i_h,:,:)), squeeze(group_av_sem(i_h,:,:)),'b');
% errorbar(1:30,squeeze(group_av_rt(i_l,:,:)), squeeze(group_av_sem(i_l,:,:)),'g'); hold off;
% xlim([0.5, 30.5]); ylim([350, 650]); xlabel('Block'); ylabel('RTs'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Overall');



% compute distractor cost
sbj_s1_drt = sbj_s1_rt([i_h,i_l],:,:) - sbj_s1_rt(i_a,:,:);
sbj_s2_drt = sbj_s2_rt([i_h,i_l],:,:) - sbj_s2_rt(i_a,:,:);
sbj_av_drt = sbj_av_rt([i_h,i_l],:,:) - sbj_av_rt(i_a,:,:);

i_h_d = 1;
i_l_d = 2;

group_s1_drt = mean(sbj_s1_drt,3);
group_s1_dsem = std(sbj_s1_drt,[],3) / sqrt(sbj_nr_minus1);
group_s2_drt = mean(sbj_s2_drt,3);
group_s2_dsem = std(sbj_s2_drt,[],3) / sqrt(sbj_nr_minus1);
group_av_drt = mean(sbj_av_drt,3);
group_av_dsem = std(sbj_av_drt,[],3) / sqrt(sbj_nr_minus1);


% plot (not in the paper)
figure; 
subplot(211); 
errorbar(1:30,squeeze(group_s1_drt(i_h_d,:)), squeeze(group_s1_dsem(i_h_d,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s1_drt(i_l_d,:)), squeeze(group_s1_dsem(i_l_d,:)),'b');
title('Session 1'); xlim([0.5, 30.5]); ylim([0, 150]); xlabel('Block'); ylabel('Distractor cost (ms)'); legend({'high','low'},'Location', 'eastoutside');
subplot(212); 
errorbar(1:30,squeeze(group_s2_drt(i_h_d,:)), squeeze(group_s2_dsem(i_h_d,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s2_drt(i_l_d,:)), squeeze(group_s2_dsem(i_l_d,:)),'b');
title('Session 2'); xlim([0.5, 30.5]); ylim([0, 150]); xlabel('Block'); ylabel('Distractor cost (ms)'); legend({'high','low'},'Location', 'eastoutside');
