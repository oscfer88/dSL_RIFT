
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


%% Get trial list

% % Condition list
% trialcond_list0 = str2num([num2str(data.target_loc) , num2str(data.distractor_loc)]); %#ok<ST2NM>
% trialcond_list0 = reshape(trialcond_list0,360,[]);
% 
% sbj_list = repmat(1:40,2,1); sbj_list = sbj_list(:); sbj_list(4) = [];
% trialcond_list = cell(40,1);
% sbj = [];
% for i=1:size(trialcond_list0,2)
%     if sbj_list(i) == sbj
%         list2 = [list; trialcond_list0(:,i)];
%         list = list2;
%         trialcond_list{sbj} = [];
%     else
%         list = trialcond_list0(:,i);
%     end
%     sbj = sbj_list(i);
%     trialcond_list{sbj} = list;
% end

% % Inter-trial repetition list
% disrepet_list0 = reshape(data.repeat_disloc,360,[]);
% sbj_list = repmat(1:40,2,1); sbj_list = sbj_list(:); sbj_list(4) = [];
% disrepet_list = cell(40,1);
% sbj = [];
% for i=1:size(disrepet_list0,2)
%     if sbj_list(i) == sbj
%         list2 = [list; disrepet_list0(:,i)];
%         list = list2;
%         disrepet_list{sbj} = [];
%     else
%         list = disrepet_list0(:,i);
%     end
%     sbj = sbj_list(i);
%     disrepet_list{sbj} = list;
% end


%% Tabulate

% [tbl,~,~,labels] = crosstab(data.subject_nr, data.distractor_loc, data.session_count) %distractor_loc
% [tbl,~,~,labels] = crosstab(data.subject_nr, data.distractor_prob, data.session_count) %distractor_prob
% [tbl,~,~,labels] = crosstab(data.distractor_loc, data.distractor_prob, data.session_count, data.subject_nr) %distractor_loc x distractor_prob

% [tbl,~,~,labels] = crosstab(data.subject_nr, data.target_loc, data.session_count) %target_loc
% [tbl,~,~,labels] = crosstab(data.subject_nr, data.target_loc, data.distractor_loc, data.session_count) %target_loc x distractor_loc


%% Remove uncorrect response trials

% data_all = data; %data = data_all;
% data = data_all(data_all.subject_nr == 20 & data_all.session_count == 2,:);
correct_trials = find(data.correct == 1);
nottoofast_trials = find(data.RT > 201);
nottooslow_trials = find(data.RT < mean(data.RT)+(2*std(data.RT)));
good_trials = intersect(correct_trials, intersect(nottoofast_trials, nottooslow_trials));
% save ('Q:\MEG_data\good_trials','good_trials')

datart = data(good_trials,:);


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
sbj_ov_rt = squeeze(mean(means2_rt,2));

sbj_s1_drt = sbj_s1_rt(2:3,:) - sbj_s1_rt(1,:);
sbj_s2_drt = sbj_s2_rt(2:3,:) - sbj_s2_rt(1,:);
sbj_ov_drt = sbj_ov_rt(2:3,:) - sbj_ov_rt(1,:);

group_s1_drt = mean(sbj_s1_drt,2);
group_s2_drt = mean(sbj_s2_drt,2);
group_ov_drt = mean(sbj_ov_drt,2);



% plot distractor cost w\ subjects
figure;

subplot(131); plot(sbj_s1_drt,'--k', 'LineWidth', .5); hold on;
plot(group_s1_drt,'k','LineWidth', 1.5); hold off;
xlim([0.5,2.5]); xlabel('Distarctor probability'); xticks(1:2); xticklabels({'high','low'});
ylim([15, 105]); ylabel('Distractor cost (ms)'); title('Day 1');
subplot(132); plot(sbj_s2_drt,'--k', 'LineWidth', .5); hold on;
plot(group_s2_drt,'k','LineWidth', 1.5); hold off;
xlim([0.5,2.5]); xlabel('Distarctor probability'); xticks(1:2); xticklabels({'high','low'});
ylim([15, 110]); ylabel('Distractor cost (ms)'); title('Day 2');
subplot(133); plot(sbj_ov_drt,'--k', 'LineWidth', .5); hold on;
plot(group_ov_drt,'k','LineWidth', 1.5); hold off;
xlim([0.5,2.5]); xlabel('Distarctor probability'); xticks(1:2); xticklabels({'high','low'});
ylim([15, 110]); ylabel('Distractor cost (ms)'); title('Overall');

set(gcf,'units','centimeters','position',[0,0,11.4,4])


% [~,p] = ttest(sbj_s1_drt(1,:), sbj_s1_drt(2,:)) %t-test
% [~,p] = ttest(sbj_s2_drt(1,:), sbj_s2_drt(2,:)) %t-test
% [~,p,~,stats] = ttest(sbj_ov_drt(1,:), sbj_ov_drt(2,:)) %t-test

% %plot posterRAW
% group_ov_dsemrt = std(sbj_ov_drt,[],2);
% group_ov_dstert = group_ov_dsemrt/sqrt(sbj_nr);
% figure; bar(group_ov_drt); hold on;
% errorbar(1:2, group_ov_drt, group_ov_dstert,'color','k','LineStyle','none');


% figure; subplot(231); plot(sbj_s1_drt,'k'); hold on; plot(group_s1_drt, 'r','LineWidth',1.5);
% xlim([.5, 2.5]); xticks([1, 2]); ylim([0, 120]); ylabel('Distractor Cost on RTs'); xticklabels({'high','low'}); title('Session 1');
% subplot(232); plot(sbj_s2_drt,'k'); hold on; plot(group_s2_drt, 'r','LineWidth',1.5);
% xlim([.5, 2.5]); xticks([1, 2]); ylim([0, 120]); xticklabels({'high','low'}); title('Session 2');
% subplot(233); plot(sbj_ov_drt,'k'); hold on; plot(group_ov_drt, 'r','LineWidth',1.5);
% xlim([.5, 2.5]); xticks([1, 2]); ylim([0, 120]); xticklabels({'high','low'}); title('Overall');
% 
% % sbj_s1_drt_m = sbj_s1_drt - sbj_s1_drt(1,:);
% % sbj_s2_drt_m = sbj_s2_drt - sbj_s2_drt(1,:);
% % sbj_ov_drt_m = sbj_ov_drt - sbj_ov_drt(1,:);
% % 
% % group_s1_drt_m = group_s1_drt - group_s1_drt(1);
% % group_s2_drt_m = group_s2_drt - group_s2_drt(1);
% % group_ov_drt_m = group_ov_drt - group_ov_drt(1);
% % 
% % figure; subplot(231); plot(sbj_s1_drt_m,'k'); hold on; plot(group_s1_drt_m, 'r','LineWidth',1.5);
% % xlim([.5, 2.5]); xticks([1, 2]); ylim([-40 40]); ylabel('Distractor Cost on RTs'); xticklabels({'high','low'}); title('Session 1');
% % subplot(232); plot(sbj_s2_drt_m,'k'); hold on; plot(group_s2_drt_m, 'r','LineWidth',1.5);
% % xlim([.5, 2.5]); xticks([1, 2]); ylim([-40 40]); xticklabels({'high','low'}); title('Session 2');
% % subplot(233); plot(sbj_ov_drt_m,'k'); hold on; plot(group_ov_drt_m, 'r','LineWidth',1.5);
% % xlim([.5, 2.5]); xticks([1, 2]); ylim([-40 40]); xticklabels({'high','low'}); title('Overall');


sbj_meanrt_ov = mean(means2_rt,2);
sbj_stdrt = mean(stds2,2);

% group_meanrt = mean(sbj_meanrt,3);
% group_stdrt = std(sbj_stdrt,[],3);
% group_stert = group_stdrt/sqrt(sbj_nr);

sbj_meanrt_s1 = means2_rt(:,1,:);
sbj_meanrt_s2 = means2_rt(:,2,:);
sbj_stdrt_s1 = stds2(:,1,:);
sbj_stdrt_s2 = stds2(:,2,:);

figure; subplot(231); boxplot(squeeze(sbj_meanrt_s1)'); hold on; plot(mean(sbj_meanrt_s1,3)','db');
ylim([375, 675]); ylabel('RTs'); xticklabels({'absent','high','low'}); title('Session 1');
subplot(232); boxplot(squeeze(sbj_meanrt_s2)'); hold on; plot(mean(sbj_meanrt_s2,3)','db');
ylim([375, 675]); xticklabels({'absent','high','low'}); title('Session 2'); %xlabel('Distractor location');
subplot(233); boxplot(squeeze(sbj_meanrt_ov)'); hold on; plot(mean(sbj_meanrt_ov,3)','db');
ylim([375, 675]); xticklabels({'absent','high','low'}); title('Overall');

% sbj_meanrt_s1_dif = squeeze(sbj_meanrt_s1(2:3,:,:)) - squeeze(sbj_meanrt_s1(1,:,:))';
% sbj_meanrt_s2_dif = squeeze(sbj_meanrt_s2(2:3,:,:)) - squeeze(sbj_meanrt_s2(1,:,:))';
% sbj_meanrt_ov_dif = squeeze(sbj_meanrt(2:3,:,:)) - squeeze(sbj_meanrt(1,:,:))';
% 
% figure;
% for i = 1:sbj_nr
%     subplot(ceil(sbj_nr/4),4,i);
%     plot(sbj_meanrt_s1_dif(:,i),'b');
%     hold on;
%     plot(sbj_meanrt_s2_dif(:,i),'g');
%     plot(sbj_meanrt_ov_dif(:,i),'r');
%     xlim([.5, 2.5]); xticks([1, 2]); ylim([20, 100]); xticklabels({'high','low'}); title(['Sbj ' num2str(i)]);
%     if ismember(i,[1 5 9 13])
%         ylabel('RTs');
%     end
%     if i == sbj_nr
%         legend('S1','S2','Ov','Location','northeast');
%     end
% end
% % figure; datasets = {sbj_meanrt_s1, sbj_meanrt_s2, sbj_meanrt}; stdevs = {sbj_stdrt_s1, sbj_stdrt_s2, sbj_stdrt};
% % for i = 1:3
% %     datasbj = datasets{i};
% %     devsbj = stdevs{i};
% %     for ii = 1:sbj_nr
% %         subplot(sbj_nr,3,i+ii+((2*ii)-3));
% %         bar([datasbj(1,1,ii),datasbj(2,1,ii),datasbj(3,1,ii)]);
% %         ylim([300, 700]); ylabel('RTs'); xticklabels({'absent','high','low'});
% %         hold on;
% %         errorbar([datasbj(1,1,ii),datasbj(2,1,ii),datasbj(3,1,ii)],[devsbj(1,1,ii),devsbj(2,1,ii),devsbj(3,1,ii)],...
% %             'LineStyle','none', 'Color', 'k');
% %     end
% % end


% group_meanrt_s1 = mean(sbj_meanrt_s1,3);
% group_meanrt_s2 = mean(sbj_meanrt_s2,3);
% group_stdrt_s1 = mean(sbj_stdrt_s1,3);
% group_stdrt_s2 = mean(sbj_stdrt_s2,3);
% group_stert_s1 = group_stdrt_s1/sqrt(sbj_nr);
% group_stert_s2 = group_stdrt_s2/sqrt(sbj_nr);
% 
% figure; subplot(231); bar([group_meanrt_s1(1),group_meanrt_s1(2),group_meanrt_s1(3)]);
% ylim([350, 650]); ylabel('RTs'); xticklabels({'absent','high','low'}); title('Session 1');
% hold on;
% errorbar([group_meanrt_s1(1),group_meanrt_s1(2),group_meanrt_s1(3)],[group_stert_s1(1),group_stert_s1(2),group_stert_s1(3)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(232); bar([group_meanrt_s2(1),group_meanrt_s2(2),group_meanrt_s2(3)]);
% ylim([350, 650]); xticklabels({'absent','high','low'}); title('Session 2'); %xlabel('Distractor location');
% hold on;
% errorbar([group_meanrt_s2(1),group_meanrt_s2(2),group_meanrt_s2(3)],[group_stert_s2(1),group_stert_s2(2),group_stert_s2(3)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(233); bar([group_meanrt(1),group_meanrt(2),group_meanrt(3)]);
% ylim([350, 650]); xticklabels({'absent','high','low'}); title('Overall');
% hold on;
% errorbar([group_meanrt(1),group_meanrt(2),group_meanrt(3)],[group_stert(1),group_stert(2),group_stert(3)],...
%     'LineStyle','none', 'Color', 'k');


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
sbj_ov_acc = squeeze(mean(means2_ac,2));

sbj_s1_dacc = sbj_s1_acc(1,:) - sbj_s1_acc(2:3,:);
sbj_s2_dacc = sbj_s2_acc(1,:) - sbj_s2_acc(2:3,:);
sbj_ov_dacc = sbj_ov_acc(1,:) - sbj_ov_acc(2:3,:);

group_s1_dacc = mean(sbj_s1_dacc,2);
group_s2_dacc = mean(sbj_s2_dacc,2);
group_ov_dacc = mean(sbj_ov_dacc,2);



% plot distractor cost w\ subjects
figure

subplot(131); plot(sbj_s1_dacc*100,'--k', 'LineWidth', .5); hold on;
plot(group_s1_dacc*100,'k','LineWidth',2); hold off;
xlim([0.5,2.5]); xlabel('Distarctor probability'); xticks(1:2); xticklabels({'high','low'});
ylim([-5, 18]); ylabel('Distractor cost (%acc)'); title('Day 1');
subplot(132); plot(sbj_s2_dacc*100,'--k'); hold on;
plot(group_s2_dacc*100,'k','LineWidth',2); hold off;
xlim([0.5,2.5]); xlabel('Distarctor probability'); xticks(1:2); xticklabels({'high','low'});
ylim([-5, 18]); ylabel('Distractor cost (%acc)'); title('Day 2');
subplot(133); plot(sbj_ov_dacc*100,'--k'); hold on;
plot(group_ov_dacc*100,'k','LineWidth',2); hold off;
xlim([0.5,2.5]); xlabel('Distarctor probability'); xticks(1:2); xticklabels({'high','low'});
ylim([-5, 18]); ylabel('Distractor cost (%acc)'); title('Overall');

set(gcf,'units','centimeters','position',[0,0,11.4,4])

% [~,p] = ttest(sbj_s1_dacc(1,:), sbj_s1_dacc(2,:)) %t-test
% [~,p] = ttest(sbj_s2_dacc(1,:), sbj_s2_dacc(2,:)) %t-test
% [~,p] = ttest(sbj_ov_dacc(1,:), sbj_ov_dacc(2,:)) %t-test

% %plot posterRAW
% group_ov_dsemacc = std(sbj_ov_dacc,[],2);
% group_ov_dsteacc = group_ov_dsemacc/sqrt(sbj_nr);
% figure; bar(group_ov_dacc); hold on;
% errorbar(1:2, group_ov_dacc, group_ov_dsteacc,'color','k','LineStyle','none');


subplot(234); plot(sbj_s1_dacc,'k'); hold on;  plot(group_s1_dacc, 'r','LineWidth',1.5);
xlim([.5, 2.5]); xticks([1, 2]); ylim([-.08, .12]); ylabel('Distractor Cost on Accuracy'); xticklabels({'high','low'}); title('Session 1');
subplot(235); plot(sbj_s2_dacc,'k'); hold on; plot(group_s2_dacc, 'r','LineWidth',1.5);
xlim([.5, 2.5]); xticks([1, 2]); ylim([-.08, .12]); xticklabels({'high','low'}); title('Session 2');
subplot(236); plot(sbj_ov_dacc,'k'); hold on; plot(group_ov_dacc, 'r','LineWidth',1.5);
xlim([.5, 2.5]); xticks([1, 2]); ylim([-.08, .12]); xticklabels({'high','low'}); title('Overall');
% 
% % sbj_s1_dacc_m = sbj_s1_dacc - sbj_s1_dacc(1,:);
% % sbj_s2_dacc_m = sbj_s2_dacc - sbj_s2_dacc(1,:);
% % sbj_ov_dacc_m = sbj_ov_dacc - sbj_ov_dacc(1,:);
% % 
% % group_s1_dacc_m = group_s1_dacc - group_s1_dacc(1);
% % group_s2_dacc_m = group_s2_dacc - group_s2_dacc(1);
% % group_ov_dacc_m = group_ov_dacc - group_ov_dacc(1);
% % 
% % subplot(234); plot(sbj_s1_dacc_m,'k'); hold on;  plot(group_s1_dacc_m, 'r','LineWidth',1.5);
% % xlim([.5, 2.5]); xticks([1, 2]); ylim([-.05, .1]); ylabel('Distractor Cost on Accuracy'); xticklabels({'high','low'}); title('Session 1');
% % subplot(235); plot(sbj_s2_dacc_m,'k'); hold on; plot(group_s2_dacc_m, 'r','LineWidth',1.5);
% % xlim([.5, 2.5]); xticks([1, 2]); ylim([-.05, .1]); xticklabels({'high','low'}); title('Session 2');
% % subplot(236); plot(sbj_ov_dacc_m,'k'); hold on; plot(group_ov_dacc_m, 'r','LineWidth',1.5);
% % xlim([.5, 2.5]); xticks([1, 2]); ylim([-.05, .1]); xticklabels({'high','low'}); title('Overall');


sbj_meanacc = mean(means2_ac,2);
sbj_stdacc = mean(stds2,2);

group_meanacc = mean(sbj_meanacc,3);
group_stdacc = std(sbj_stdacc,[],3);
group_steacc = group_stdacc/sqrt(sbj_nr);

sbj_meanacc_s1 = means2_ac(:,1,:);
sbj_meanacc_s2 = means2_ac(:,2,:);
sbj_stdacc_s1 = stds2(:,1,:);
sbj_stdacc_s2 = stds2(:,2,:);

subplot(234); boxplot(squeeze(sbj_meanacc_s1([1 2 3],:,:))'); hold on; plot(mean(sbj_meanacc_s1([2 1 3],:,:),3)','db');
ylim([.8, 1]); ylabel('Accuracy'); xticklabels({'absent','high','low'}); title('Session 1');
subplot(235); boxplot(squeeze(sbj_meanacc_s2([1 2 3],:,:))'); hold on; plot(mean(sbj_meanacc_s2([2 1 3],:,:),3)','db');
ylim([.8, 1]); xticklabels({'absent','high','low'}); title('Session 2'); %xlabel('Distractor location');
subplot(236); boxplot(squeeze(sbj_meanacc([1 2 3],:,:))'); hold on; plot(mean(sbj_meanacc([2 1 3],:,:),3)','db');
ylim([.8, 1]); xticklabels({'absent','high','low'}); title('Overall');


% sbj_meanacc_s1_dif = squeeze(sbj_meanacc_s1(2,:,:))' - squeeze(sbj_meanacc_s1([1 3],:,:));
% sbj_meanacc_s2_dif = squeeze(sbj_meanacc_s2(2,:,:))' - squeeze(sbj_meanacc_s2([1 3],:,:));
% sbj_meanacc_ov_dif = squeeze(sbj_meanacc(2,:,:))' - squeeze(sbj_meanacc([1 3],:,:));
% 
% figure;
% for i = 1:sbj_nr
%     subplot(ceil(sbj_nr/4),4,i);
%     plot(sbj_meanacc_s1_dif(:,i),'b');
%     hold on;
%     plot(sbj_meanacc_s2_dif(:,i),'g');
%     plot(sbj_meanacc_ov_dif(:,i),'r');
%     xlim([.5, 2.5]); xticks([1, 2]); ylim([-.04 .16]); xticklabels({'high','low'});
%     if ismember(i,[1 5 9 13])
%         ylabel('Accuracy');
%     end
%     if i == sbj_nr
%         legend('S1','S2','Ov','Location','northeast');
%     end
% end
% % figure; datasets = {sbj_meanacc_s1, sbj_meanacc_s2, sbj_meanacc}; stdevs = {sbj_stdacc_s1, sbj_stdacc_s2, sbj_stdacc};
% % for i = 1:3
% %     datasbj = datasets{i};
% %     devsbj = stdevs{i};
% %     for ii = 1:4
% %         subplot(4,3,i+ii+((2*ii)-3));
% %         bar([datasbj(2,1,ii),datasbj(1,1,ii),datasbj(3,1,ii)]);
% %         ylim([.7, 1]); ylabel('Accuracy'); xticklabels({'absent','high','low'});
% %         hold on;
% %         errorbar([datasbj(1,1,ii),datasbj(2,1,ii),datasbj(3,1,ii)],[devsbj(1,1,ii),devsbj(2,1,ii),devsbj(3,1,ii)],...
% %             'LineStyle','none', 'Color', 'k');
% %     end
% % end


% group_meanacc_s1 = mean(sbj_meanacc_s1,3);
% group_meanacc_s2 = mean(sbj_meanacc_s2,3);
% group_stdacc_s1 = mean(sbj_stdacc_s1,3);
% group_stdacc_s2 = mean(sbj_stdacc_s2,3);
% group_steacc_s1 = group_stdacc_s1/sqrt(sbj_nr);
% group_steacc_s2 = group_stdacc_s2/sqrt(sbj_nr);
% 
% subplot(234); bar([group_meanacc_s1(2),group_meanacc_s1(1),group_meanacc_s1(3)]);
% ylim([.8, 1]); ylabel('Accuracy'); xticklabels({'absent','high','low'}); %title('Session 1');
% hold on;
% errorbar([group_meanacc_s1(2),group_meanacc_s1(1),group_meanacc_s1(3)],[group_steacc_s1(2),group_steacc_s1(1),group_steacc_s1(3)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(235); bar([group_meanacc_s2(2),group_meanacc_s2(1),group_meanacc_s2(3)]);
% ylim([.8, 1]); xlabel('Distractor location'); xticklabels({'absent','high','low'}); %title('Session 2');
% hold on;
% errorbar([group_meanacc_s2(2),group_meanacc_s2(1),group_meanacc_s2(3)],[group_steacc_s2(2),group_steacc_s2(1),group_steacc_s2(3)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(236); bar([group_meanacc(2),group_meanacc(1),group_meanacc(3)]);
% ylim([.8, 1]); xticklabels({'absent','high','low'}); %title('Overall');
% hold on;
% errorbar([group_meanacc(2),group_meanacc(1),group_meanacc(3)],[group_steacc(2),group_steacc(1),group_steacc(3)],...
%     'LineStyle','none', 'Color', 'k');


%% Distractor location - IES analysis

sbj_s1_ies = sbj_s1_rt ./ sbj_s1_acc;
sbj_s2_ies = sbj_s2_rt ./ sbj_s2_acc;
% sbj_ov_ies = sbj_ov_rt ./ sbj_ov_acc;
sbj_ov_ies = mean(cat(3,sbj_s1_ies, sbj_s2_ies),3);

sbj_s1_dies = sbj_s1_ies(2:3,:) - sbj_s1_ies(1,:);
sbj_s2_dies = sbj_s2_ies(2:3,:) - sbj_s2_ies(1,:);
sbj_ov_dies = sbj_ov_ies(2:3,:) - sbj_ov_ies(1,:);


% % individual plot
% x = 10;
% figure; subplot(131); bar(sbj_s1_ies(:,x)); ylim([350, 750]); ylabel('IES'); xticklabels({'absent','high','low'}); title('Session 1');
% subplot(132); bar(sbj_s2_ies(:,x)); ylim([350, 750]); xlabel('Distractor location');xticklabels({'absent','high','low'}); title('Session 2');
% subplot(133); bar(sbj_ov_ies(:,x)); ylim([350, 750]); xticklabels({'absent','high','low'}); title('Overall');


group_s1_dies = mean(sbj_s1_dies,2);
group_s2_dies = mean(sbj_s2_dies,2);
group_ov_dies = mean(sbj_ov_dies,2);

% [~,p] = ttest(sbj_s1_dies(1,:), sbj_s1_dies(2,:)) %t-test
% [~,p] = ttest(sbj_s2_dies(1,:), sbj_s2_dies(2,:)) %t-test
% [~,p] = ttest(sbj_ov_dies(1,:), sbj_ov_dies(2,:)) %t-test

% figure; subplot(131); plot(sbj_s1_dies,'k'); hold on;  plot(group_s1_dies, 'r','LineWidth',1.5);
% xlim([.5, 2.5]); xticks([1, 2]); ylim([0 150]); ylabel('Distractor Cost on Acc (ms)'); xticklabels({'high','low'}); title('Session 1');
% subplot(132); plot(sbj_s2_dies,'k'); hold on; plot(group_s2_dies, 'r','LineWidth',1.5);
% xlim([.5, 2.5]); xticks([1, 2]); ylim([0 150]); xticklabels({'high','low'}); title('Session 2');
% subplot(133); plot(sbj_ov_dies,'k'); hold on; plot(group_ov_dies, 'r','LineWidth',1.5);
% xlim([.5, 2.5]); xticks([1, 2]); ylim([0 150]); xticklabels({'high','low'}); title('Overall');
% 
% % sbj_s1_dies_m = sbj_s1_dies - sbj_s1_dies(1,:);
% % sbj_s2_dies_m = sbj_s2_dies - sbj_s2_dies(1,:);
% % sbj_ov_dies_m = sbj_ov_dies - sbj_ov_dies(1,:);
% % 
% % group_s1_dies_m = group_s1_dies - group_s1_dies(1);
% % group_s2_dies_m = group_s2_dies - group_s2_dies(1);
% % group_ov_dies_m = group_ov_dies - group_ov_dies(1);
% % 
% % figure; subplot(131); plot(sbj_s1_dies_m,'k'); hold on;  plot(group_s1_dies_m, 'r','LineWidth',1.5);
% % xlim([.5, 2.5]); xticks([1, 2]); ylim([-50, 130]); ylabel('Distractor Cost on IES'); xticklabels({'high','low'}); title('Session 1');
% % subplot(132); plot(sbj_s2_dies_m,'k'); hold on; plot(group_s2_dies_m, 'r','LineWidth',1.5);
% % xlim([.5, 2.5]); xticks([1, 2]); ylim([-50, 130]); xticklabels({'high','low'}); title('Session 2');
% % subplot(133); plot(sbj_ov_dies_m,'k'); hold on; plot(group_ov_dies_m, 'r','LineWidth',1.5);
% % xlim([.5, 2.5]); xticks([1, 2]); ylim([-50, 130]); xticklabels({'high','low'}); title('Overall');


sbj_meanies = sbj_meanrt_ov ./ sbj_meanacc;
sbj_meanies_s1 = sbj_meanrt_s1 ./ sbj_meanacc_s1;
sbj_meanies_s2 = sbj_meanrt_s2 ./ sbj_meanacc_s2;

figure; subplot(131); boxplot(squeeze(sbj_meanies_s1)'); hold on; plot(mean(sbj_meanies_s1,3)','db');
ylim([400, 700]); ylabel('IES'); xticklabels({'absent','high','low'}); title('Session 1');
subplot(132); boxplot(squeeze(sbj_meanies_s2)'); hold on; plot(mean(sbj_meanies_s2,3)','db');
ylim([400, 700]); xticklabels({'absent','high','low'}); title('Session 2'); %xlabel('Distractor location');
subplot(133); boxplot(squeeze(sbj_meanies)'); hold on; plot(mean(sbj_meanies,3)','db');
ylim([400, 700]); xticklabels({'absent','high','low'}); title('Overall');


% sbj_meanies_s1_dif = squeeze(sbj_meanies_s1(2:3,:,:)) - squeeze(sbj_meanies_s1(1,:,:))';
% sbj_meanies_s2_dif = squeeze(sbj_meanies_s2(2:3,:,:)) - squeeze(sbj_meanies_s2(1,:,:))';
% sbj_meanies_ov_dif = squeeze(sbj_meanies(2:3,:,:)) - squeeze(sbj_meanies(1,:,:))';
% 
% figure;
% for i = 1:sbj_nr
%     subplot(ceil(sbj_nr/4),4,i);
%     plot(sbj_meanies_s1_dif(:,i),'b');
%     hold on;
%     plot(sbj_meanies_s2_dif(:,i),'g');
%     plot(sbj_meanies_ov_dif(:,i),'r');
%     xlim([.5, 2.5]); xticks([1, 2]); ylim([0 170]); xticklabels({'high','low'});
%     if ismember(i,[1 5 9 13])
%         ylabel('IES');
%     end
%     if i == sbj_nr
%         legend('S1','S2','Ov','Location','northeast');
%     end
% end


% group_meanies_s1 = group_meanrt_s1 ./ group_meanacc_s1;
% group_meanies_s2 = group_meanrt_s2 ./ group_meanacc_s2;
% group_meanies = group_meanrt ./ group_meanacc;
% 
% figure; subplot(131); bar([group_meanies_s1(1),group_meanies_s1(2),group_meanies_s1(3)]);
% ylim([350, 650]); ylabel('IES'); xticklabels({'absent','high','low'}); title('Session 1');
% subplot(132); bar([group_meanies_s2(1),group_meanies_s2(2),group_meanies_s2(3)]);
% ylim([350, 650]); xticklabels({'absent','high','low'}); title('Session 2'); %xlabel('Distractor location');
% subplot(133); bar([group_meanies(1),group_meanies(2),group_meanies(3)]);
% ylim([350, 650]); xticklabels({'absent','high','low'}); title('Overall');


%% Distractor location - Balanced integration score (BIS) analysis          DOULBE-CHECK THIS!!!!!!!!!!!!!!!!!

sbj_s1_zrt = zscore(sbj_s1_rt);
sbj_s1_zacc = zscore(sbj_s1_acc);
sbj_s2_zrt = zscore(sbj_s2_rt);
sbj_s2_zacc = zscore(sbj_s2_acc);

sbj_s1_bis = sbj_s1_zacc - sbj_s1_zrt;
sbj_s2_bis = sbj_s2_zacc - sbj_s2_zrt;
sbj_ov_bis = mean(cat(3,sbj_s1_bis, sbj_s2_bis),3);

% [~,p] = ttest(sbj_s1_bis(2,:), sbj_s1_bis(3,:)) %t-test
% [~,p] = ttest(sbj_s2_bis(2,:), sbj_s2_bis(3,:)) %t-test
% [~,p] = ttest(sbj_ov_bis(2,:), sbj_ov_bis(3,:)) %t-test

figure; subplot(131); boxplot(sbj_s1_bis'); hold on; plot(mean(sbj_s1_bis,2)','db');
ylim([-2.5, 2.5]); ylabel('BIS'); xticklabels({'absent','high','low'}); title('Session 1');
subplot(132); boxplot(sbj_s2_bis'); hold on; plot(mean(sbj_s2_bis,2)','db');
ylim([-2.5, 2.5]); xticklabels({'absent','high','low'}); title('Session 2'); %xlabel('Distractor location');
subplot(133); boxplot(sbj_ov_bis'); hold on; plot(mean(sbj_ov_bis,2)','db');
ylim([-2.5, 2.5]); xticklabels({'absent','high','low'}); title('Overall');


%% Distractor location x Block- RT analysis
% Mean reaction time
[means,stds,grp] = grpstats(datart.RT,...
    {datart.subject_nr, datart.session_count, datart.block_count, datart.distractor_prob},...
    {'mean','std','gname'});
means2 = reshape(means,3,30,2,sbj_nr);
stds2 = reshape(stds,3,30,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),'_',grp(:,2),'_',grp(:,3),'_',grp(:,4)),3,30,2,sbj_nr);

sbj_s1_rt = squeeze(means2(:,:,1,:));
sbj_s2_rt = squeeze(means2(:,:,2,:));
sbj_ov_rt = squeeze(mean(means2,3));

group_s1_rt = mean(sbj_s1_rt,3);
group_s1_sem = std(sbj_s1_rt,[],3) / sqrt(sbj_nr);
group_s2_rt = mean(sbj_s2_rt,3);
group_s2_sem = std(sbj_s2_rt,[],3) / sqrt(sbj_nr);
group_ov_rt = mean(sbj_ov_rt,3);
group_ov_sem = std(sbj_ov_rt,[],3) / sqrt(sbj_nr);

figure;
subplot(311);
errorbar(1:30,squeeze(group_s1_rt(1,:,:)), squeeze(group_s1_sem(1,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s1_rt(2,:,:)), squeeze(group_s1_sem(2,:,:)),'b');
errorbar(1:30,squeeze(group_s1_rt(3,:,:)), squeeze(group_s1_sem(3,:,:)),'g');
xlim([0.5, 30.5]); ylim([350, 650]); ylabel('RTs'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Session 1');
subplot(312);
errorbar(1:30,squeeze(group_s2_rt(1,:,:)), squeeze(group_s2_sem(1,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s2_rt(2,:,:)), squeeze(group_s2_sem(2,:,:)),'b');
errorbar(1:30,squeeze(group_s2_rt(3,:,:)), squeeze(group_s2_sem(3,:,:)),'g');
xlim([0.5, 30.5]); ylim([350, 650]); ylabel('RTs'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Session 2');
subplot(313);
errorbar(1:30,squeeze(group_ov_rt(1,:,:)), squeeze(group_ov_sem(1,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_ov_rt(2,:,:)), squeeze(group_ov_sem(2,:,:)),'b');
errorbar(1:30,squeeze(group_ov_rt(3,:,:)), squeeze(group_ov_sem(3,:,:)),'g'); hold off;
xlim([0.5, 30.5]); ylim([430, 560]); xlabel('Block'); ylabel('RTs'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Overall');



%distractor cost
sbj_s1_drt = sbj_s1_rt(2:3,:,:) - sbj_s1_rt(1,:,:);
sbj_s2_drt = sbj_s2_rt(2:3,:,:) - sbj_s2_rt(1,:,:);
sbj_ov_drt = sbj_ov_rt(2:3,:,:) - sbj_ov_rt(1,:,:);

group_s1_drt = mean(sbj_s1_drt,3);
group_s1_dsem = std(sbj_s1_drt,[],3) / sqrt(sbj_nr);
group_s2_drt = mean(sbj_s2_drt,3);
group_s2_dsem = std(sbj_s2_drt,[],3) / sqrt(sbj_nr);
group_ov_drt = mean(sbj_ov_drt,3);
group_ov_dsem = std(sbj_ov_drt,[],3) / sqrt(sbj_nr);

figure; subplot(311); errorbar(1:30,squeeze(group_s1_drt(1,:)), squeeze(group_s1_dsem(1,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s1_drt(2,:)), squeeze(group_s1_dsem(2,:)),'b');
title('Day 1'); xlim([0.5, 30.5]); ylim([0, 110]); xlabel('Block'); %ylabel('Distractor cost (ms)'); %legend({'high','low'},'Location', 'eastoutside');
subplot(312); errorbar(1:30,squeeze(group_s2_drt(1,:)), squeeze(group_s2_dsem(1,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s2_drt(2,:)), squeeze(group_s2_dsem(2,:)),'b');
title('Day 2'); xlim([0.5, 30.5]); ylim([0, 110]); xlabel('Block'); ylabel('Distractor cost (ms)'); %legend({'high','low'},'Location', 'eastoutside');
subplot(313); errorbar(1:30,squeeze(group_ov_drt(1,:)), squeeze(group_ov_dsem(1,:)),'r'); hold on;
errorbar(1:30,squeeze(group_ov_drt(2,:)), squeeze(group_ov_dsem(2,:)),'b');
title('Overall'); xlim([0.5, 30.5]); ylim([0, 110]); xlabel('Block'); %ylabel('Distractor cost (ms)'); %legend({'high','low'},'Location', 'eastoutside');






%% Distractor location x Block- Accuracy analysis
% Mean accuracy
[means,stds,grp] = grpstats(data.correct,...
    {data.subject_nr, data.session_count, data.block_count, data.distractor_prob},...
    {'mean','std','gname'});
means2 = reshape(means,3,30,2,sbj_nr);
stds2 = reshape(stds,3,30,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),'_',grp(:,2),'_',grp(:,3),'_',grp(:,4)),3,30,2,sbj_nr);

sbj_s1_acc = squeeze(means2(:,:,1,:));
sbj_s2_acc = squeeze(means2(:,:,2,:));
sbj_ov_acc = squeeze(mean(means2,3));

% sbj_s1_dacc = sbj_s1_acc(2,:,:) - sbj_s1_acc([1 3],:,:);
% sbj_s2_dacc = sbj_s2_acc(2,:,:) - sbj_s2_acc([1 3],:,:);
% sbj_ov_dacc = sbj_ov_acc(2,:,:) - sbj_ov_acc([1 3],:,:);

group_s1_acc = mean(sbj_s1_acc,3);
group_s1_sem = std(sbj_s1_acc,[],3) / sqrt(sbj_nr);
group_s2_acc = mean(sbj_s2_acc,3);
group_s2_sem = std(sbj_s2_acc,[],3) / sqrt(sbj_nr);
group_ov_acc = mean(sbj_ov_acc,3);
group_ov_sem = std(sbj_ov_acc,[],3) / sqrt(sbj_nr);

figure;
subplot(311);
errorbar(1:30,squeeze(group_s1_acc(2,:,:)), squeeze(group_s1_sem(2,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s1_acc(1,:,:)), squeeze(group_s1_sem(1,:,:)),'b');
errorbar(1:30,squeeze(group_s1_acc(3,:,:)), squeeze(group_s1_sem(3,:,:)),'g');
xlim([0.5, 30.5]); ylim([0.8 1]); ylabel('Accuracy'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Session 1');
subplot(312);
errorbar(1:30,squeeze(group_s2_acc(2,:,:)), squeeze(group_s2_sem(2,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s2_acc(1,:,:)), squeeze(group_s2_sem(1,:,:)),'b');
errorbar(1:30,squeeze(group_s2_acc(3,:,:)), squeeze(group_s2_sem(3,:,:)),'g');
xlim([0.5, 30.5]); ylim([0.8 1]); ylabel('Accuracy'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Session 2');
subplot(313);
errorbar(1:30,squeeze(group_ov_acc(2,:,:)), squeeze(group_ov_sem(2,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_ov_acc(1,:,:)), squeeze(group_ov_sem(1,:,:)),'b');
errorbar(1:30,squeeze(group_ov_acc(3,:,:)), squeeze(group_ov_sem(3,:,:)),'g');
xlim([0.5, 30.5]); ylim([0.8 1]); ylabel('Accuracy'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Overall');


%% Distractor location x Block- IES analysis

sbj_s1_ies = sbj_s1_rt ./ sbj_s1_acc([2 1 3],:,:);
sbj_s2_ies = sbj_s2_rt ./ sbj_s2_acc([2 1 3],:,:);
sbj_ov_ies = sbj_ov_rt ./ sbj_ov_acc([2 1 3],:,:);

group_s1_ies = mean(sbj_s1_ies,3);
group_s1_sem = std(sbj_s1_ies,[],3) / sqrt(sbj_nr);
group_s2_ies = mean(sbj_s2_ies,3);
group_s2_sem = std(sbj_s2_ies,[],3) / sqrt(sbj_nr);
group_ov_ies = mean(sbj_ov_ies,3);
group_ov_sem = std(sbj_ov_ies,[],3) / sqrt(sbj_nr);

figure;
subplot(311);
errorbar(1:30,squeeze(group_s1_ies(1,:,:)), squeeze(group_s1_sem(1,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s1_ies(2,:,:)), squeeze(group_s1_sem(2,:,:)),'b');
errorbar(1:30,squeeze(group_s1_ies(3,:,:)), squeeze(group_s1_sem(3,:,:)),'g');
xlim([0.5, 30.5]); ylim([400, 700]); ylabel('IES'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Session 1');
subplot(312);
errorbar(1:30,squeeze(group_s2_ies(1,:,:)), squeeze(group_s2_sem(1,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_s2_ies(2,:,:)), squeeze(group_s2_sem(2,:,:)),'b');
errorbar(1:30,squeeze(group_s2_ies(3,:,:)), squeeze(group_s2_sem(3,:,:)),'g');
xlim([0.5, 30.5]); ylim([400, 700]); ylabel('IES'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Session 2');
subplot(313);
errorbar(1:30,squeeze(group_ov_ies(1,:,:)), squeeze(group_ov_sem(1,:,:)),'r'); hold on;
errorbar(1:30,squeeze(group_ov_ies(2,:,:)), squeeze(group_ov_sem(2,:,:)),'b');
errorbar(1:30,squeeze(group_ov_ies(3,:,:)), squeeze(group_ov_sem(3,:,:)),'g');
xlim([0.5, 30.5]); ylim([400, 700]); ylabel('IES'); legend({'absent','high','low'},'Location', 'eastoutside'); title('Overall');


%% Inter-trial distractor location repetition - RTs analysis

% ditractor location effect without repetitions
[means,stds,grp] = grpstats(datart.RT(datart.repeat_disloc == 'n'),...
    {datart.subject_nr(datart.repeat_disloc == 'n'), datart.session_count(datart.repeat_disloc == 'n'), datart.distractor_prob(datart.repeat_disloc == 'n')},...
    {'mean','std','gname'});
means2_rt = reshape(means,3,2,sbj_nr);
stds2 = reshape(stds,3,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),3,2,sbj_nr);

sbj_s1_rt = squeeze(means2_rt(:,1,:));
sbj_s2_rt = squeeze(means2_rt(:,2,:));
sbj_ov_rt = squeeze(mean(means2_rt,2));

sbj_s1_drt = sbj_s1_rt(2:3,:) - sbj_s1_rt(1,:);
sbj_s2_drt = sbj_s2_rt(2:3,:) - sbj_s2_rt(1,:);
sbj_ov_drt = sbj_ov_rt(2:3,:) - sbj_ov_rt(1,:);

% [~,p] = ttest(sbj_s1_drt(1,:), sbj_s1_drt(2,:)) %t-test
% [~,p] = ttest(sbj_s2_drt(1,:), sbj_s2_drt(2,:)) %t-test
% [~,p] = ttest(sbj_ov_drt(1,:), sbj_ov_drt(2,:)) %t-test

group_s1_drt = mean(sbj_s1_drt,2);
group_s2_drt = mean(sbj_s2_drt,2);
group_ov_drt = mean(sbj_ov_drt,2);

figure; subplot(231); boxplot(sbj_s1_drt'); hold on; plot(group_s1_drt','db');
ylim([0, 120]); ylabel('RTs'); xticklabels({'high','low'}); title('Session 1'); xlabel('Distractor location');
subplot(232); boxplot(sbj_s2_drt'); hold on; plot(group_s2_drt','db');
ylim([0, 120]); xticklabels({'high','low'}); title('Session 2'); xlabel('Distractor location');
subplot(233); boxplot(sbj_ov_drt'); hold on; plot(group_ov_drt','db');
ylim([0, 120]); xticklabels({'high','low'}); title('Overall'); xlabel('Distractor location');

% inter-trial priming effect
[means,stds,grp] = grpstats(datart.RT,...
    {datart.subject_nr, datart.session_count, datart.repeat_disloc},...
    {'mean','std','gname'});
means2_rt = reshape(means,2,2,sbj_nr);
stds2 = reshape(stds,2,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),2,2,sbj_nr);

sbj_s1_rt = squeeze(means2_rt(:,1,:));
sbj_s2_rt = squeeze(means2_rt(:,2,:));
sbj_ov_rt = squeeze(mean(means2_rt,2));

% [~,p] = ttest(sbj_s1_rt(1,:), sbj_s1_rt(2,:)) %t-test
% [~,p] = ttest(sbj_s2_rt(1,:), sbj_s2_rt(2,:)) %t-test
% [~,p] = ttest(sbj_ov_rt(1,:), sbj_ov_rt(2,:)) %t-test

group_s1_rt = mean(sbj_s1_rt,2);
group_s2_rt = mean(sbj_s2_rt,2);
group_ov_rt = mean(sbj_ov_rt,2);

subplot(234); boxplot(sbj_s1_rt'); hold on; plot(group_s1_rt','db');
ylim([350, 650]); ylabel('RTs'); xticklabels({'no','yes'}); title('Session 1 - Inter-trial priming'); xlabel('Distractor loc repetition');
subplot(235); boxplot(sbj_s2_rt'); hold on; plot(group_s2_rt','db');
ylim([350, 650]); xticklabels({'no','yes'}); title('Session 2 - Inter-trial priming'); xlabel('Distractor loc repetition');
subplot(236); boxplot(sbj_ov_rt'); hold on; plot(group_ov_rt','db');
ylim([350, 650]); xticklabels({'no','yes'}); title('Overall - Inter-trial priming'); xlabel('Distractor loc repetition');


%% Target location - RT analysis
% Mean reaction time
[means,stds,grp] = grpstats(datart.RT,...
    {datart.subject_nr, datart.session_count, datart.target_prob},...
    {'mean','std','gname'});
means2 = reshape(means,2,2,sbj_nr);
stds2 = reshape(stds,2,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),2,2,sbj_nr);

sbj_meanrt_s1 = squeeze(means2(:,1,:));
sbj_meanrt_s2 = squeeze(means2(:,2,:));
sbj_meanrt_ov = squeeze(mean(means2,2));
% sbj_stdrt_s1 = squeeze(stds2(:,1,:));
% sbj_stdrt_s2 = squeeze(stds2(:,2,:));

% [~,p] = ttest(sbj_meanrt_s1(1,:), sbj_meanrt_s1(2,:)) %t-test
% [~,p] = ttest(sbj_meanrt_s2(1,:), sbj_meanrt_s2(2,:)) %t-test
% [~,p] = ttest(sbj_meanrt_ov(1,:), sbj_meanrt_ov(2,:)) %t-test

figure; subplot(231); boxplot(sbj_meanrt_s1'); hold on; plot(mean(sbj_meanrt_s1,2)','db');
ylim([350, 650]); ylabel('RTs'); xticklabels({'high','low'}); title('Session 1');
subplot(232); boxplot(sbj_meanrt_s2'); hold on; plot(mean(sbj_meanrt_s2,2)','db');
ylim([350, 650]); xticklabels({'high','low'}); title('Session 2'); %xlabel('Target location');
subplot(233); boxplot(sbj_meanrt_ov'); hold on; plot(mean(sbj_meanrt_ov,2)','db');
ylim([350, 650]); xticklabels({'high','low'}); title('Overall');

% figure; datasets = {sbj_meanrt_s1, sbj_meanrt_s2, sbj_meanrt}; stdevs = {sbj_stdrt_s1, sbj_stdrt_s2, sbj_stdrt};
% for i = 1:3
%     datasbj = datasets{i};
%     devsbj = stdevs{i};
%     for ii = 1:sbj_nr
%         subplot(sbj_nr,3,i+ii+((2*ii)-3));
%         bar([datasbj(1,1,ii),datasbj(2,1,ii)]);
%         ylim([300, 700]); ylabel('RTs'); xticklabels({'high','low'});
%         hold on;
%         errorbar([datasbj(1,1,ii),datasbj(2,1,ii)],[devsbj(1,1,ii),devsbj(2,1,ii)],...
%             'LineStyle','none', 'Color', 'k');
%     end
% end

% group_meanrt_s1 = mean(sbj_meanrt_s1,3);
% group_meanrt_s2 = mean(sbj_meanrt_s2,3);
% group_stdrt_s1 = mean(sbj_stdrt_s1,3);
% group_stdrt_s2 = mean(sbj_stdrt_s2,3);
% group_stert_s1 = group_stdrt_s1/sqrt(sbj_nr);
% group_stert_s2 = group_stdrt_s2/sqrt(sbj_nr);
% 
% figure; subplot(231); bar([group_meanrt_s1(1),group_meanrt_s1(2)]);
% ylim([350, 650]); ylabel('RTs'); xticklabels({'high','low'}); title('Session 1');
% hold on;
% errorbar([group_meanrt_s1(1),group_meanrt_s1(2)],[group_stert_s1(1),group_stert_s1(2)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(232); bar([group_meanrt_s2(1),group_meanrt_s2(2)]);
% ylim([350, 650]); xticklabels({'high','low'}); title('Session 2'); %xlabel('Target location'); 
% hold on;
% errorbar([group_meanrt_s2(1),group_meanrt_s2(2)],[group_stert_s2(1),group_stert_s2(2)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(233); bar([group_meanrt(1),group_meanrt(2)]);
% ylim([350, 650]); xticklabels({'high','low'}); title('Overall');
% hold on;
% errorbar([group_meanrt(1),group_meanrt(2)],[group_stert(1),group_stert(2)],...
%     'LineStyle','none', 'Color', 'k');


%% Target location - Accuracy analysis
% Mean accuracy
[means,stds,grp] = grpstats(data.correct,...
    {data.subject_nr, data.session_count, data.target_prob},...
    {'mean','std','gname'});
means2 = reshape(means,2,2,sbj_nr);
stds2 = reshape(stds,2,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),2,2,sbj_nr);

sbj_meanacc_s1 = squeeze(means2(:,1,:));
sbj_meanacc_s2 = squeeze(means2(:,2,:));
sbj_meanacc_ov = squeeze(mean(means2,2));
% sbj_stdacc_s1 = stds2(:,1,:);
% sbj_stdacc_s2 = stds2(:,2,:);

% [~,p] = ttest(sbj_meanacc_s1(1,:), sbj_meanacc_s1(2,:)) %t-test
% [~,p] = ttest(sbj_meanacc_s2(1,:), sbj_meanacc_s2(2,:)) %t-test
% [~,p] = ttest(sbj_meanacc_ov(1,:), sbj_meanacc_ov(2,:)) %t-test

subplot(234); boxplot(sbj_meanacc_s1'); hold on; plot(mean(sbj_meanacc_s1,2)','db');
ylim([.8, 1]); ylabel('Accuracy'); xticklabels({'high','low'}); title('Session 1');
subplot(235); boxplot(sbj_meanacc_s2'); hold on; plot(mean(sbj_meanacc_s2,2)','db');
ylim([.8, 1]); xticklabels({'high','low'}); title('Session 2'); %xlabel('Target location');
subplot(236); boxplot(sbj_meanacc_ov'); hold on; plot(mean(sbj_meanacc_ov,2)','db');
ylim([.8, 1]); xticklabels({'high','low'}); title('Overall');

% group_meanacc_s1 = mean(sbj_meanacc_s1,3);
% group_meanacc_s2 = mean(sbj_meanacc_s2,3);
% group_stdacc_s1 = mean(sbj_stdacc_s1,3);
% group_stdacc_s2 = mean(sbj_stdacc_s2,3);
% group_steacc_s1 = group_stdacc_s1/sqrt(sbj_nr);
% group_steacc_s2 = group_stdacc_s2/sqrt(sbj_nr);
% 
% subplot(234); bar([group_meanacc_s1(1),group_meanacc_s1(2)]);
% ylim([.8, 1]); ylabel('Accuracy'); xticklabels({'high','low'}); %title('Session 1');
% hold on;
% errorbar([group_meanacc_s1(1),group_meanacc_s1(2)],[group_steacc_s1(1),group_steacc_s1(2)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(235); bar([group_meanacc_s2(1),group_meanacc_s2(2)]);
% ylim([.8, 1]); xlabel('Target location'); xticklabels({'high','low'}); %title('Session 2');
% hold on;
% errorbar([group_meanacc_s2(1),group_meanacc_s2(2)],[group_steacc_s2(1),group_steacc_s2(2)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(236); bar([group_meanacc(1),group_meanacc(2)]);
% ylim([.8, 1]); xticklabels({'high','low'}); %title('Overall');
% hold on;
% errorbar([group_meanacc(1),group_meanacc(2)],[group_steacc(1),group_steacc(2)],...
%     'LineStyle','none', 'Color', 'k');


%% Target location - IES analysis

sbj_meanies_s1 = sbj_meanrt_s1 ./ sbj_meanacc_s1;
sbj_meanies_s2 = sbj_meanrt_s2 ./ sbj_meanacc_s2;
% sbj_meanies_ov = sbj_meanrt_ov ./ sbj_meanacc_ov;
sbj_meanies_ov = mean(cat(3,sbj_meanies_s1, sbj_meanies_s2),3);

% [~,p] = ttest(sbj_meanies_s1(1,:), sbj_meanies_s1(2,:)) %t-test
% [~,p] = ttest(sbj_meanies_s2(1,:), sbj_meanies_s2(2,:)) %t-test
% [~,p] = ttest(sbj_meanies_ov(1,:), sbj_meanies_ov(2,:)) %t-test

figure; subplot(131); boxplot(sbj_meanies_s1'); hold on; plot(mean(sbj_meanies_s1,2)','db');
ylim([400, 700]); ylabel('IES'); xticklabels({'high','low'}); title('Session 1');
subplot(132); boxplot(sbj_meanies_s2'); hold on; plot(mean(sbj_meanies_s2,2)','db');
ylim([400, 700]); xticklabels({'high','low'}); title('Session 2'); %xlabel('Target location');
subplot(133); boxplot(sbj_meanies_ov'); hold on; plot(mean(sbj_meanies_ov,2)','db');
ylim([400, 700]); xticklabels({'high','low'}); title('Overall');

% group_meanies_s1 = group_meanrt_s1 ./ group_meanacc_s1;
% group_meanies_s2 = group_meanrt_s2 ./ group_meanacc_s2;
% group_meanies = group_meanrt ./ group_meanacc;
% 
% figure; subplot(131); bar([group_meanies_s1(1),group_meanies_s1(2)]);
% ylim([350, 650]); ylabel('IES'); xticklabels({'high','low'}); title('Session 1');
% subplot(132); bar([group_meanies_s2(1),group_meanies_s2(2)]);
% ylim([350, 650]); xticklabels({'high','low'}); title('Session 2'); %xlabel('Distractor location');
% subplot(133); bar([group_meanies(1),group_meanies(2)]);
% ylim([350, 650]); xticklabels({'high','low'}); title('Overall');


%% Target location X Distractor Presence - RT analysis
% Mean reaction time
[means,stds,grp] = grpstats(datart.RT,...
    {datart.subject_nr, datart.session_count, datart.target_prob, datart.distractor_pres},...
    {'mean','std','gname'});
means2 = reshape(means,2,2,2,sbj_nr);
stds2 = reshape(stds,2,2,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3),grp(:,4)),2,2,2,sbj_nr);

sbj_meanrt_s1_da = squeeze(means2(1,:,1,:)); %distractor absent
sbj_meanrt_s2_da = squeeze(means2(1,:,2,:));
sbj_meanrt_ov_da = squeeze(mean(means2(1,:,:,:),3));
sbj_meanrt_s1_dp = squeeze(means2(2,:,1,:)); %distractor present
sbj_meanrt_s2_dp = squeeze(means2(2,:,2,:));
sbj_meanrt_ov_dp = squeeze(mean(means2(2,:,:,:),3));

% [~,p] = ttest(sbj_meanrt_s1_da(1,:), sbj_meanrt_s1_da(2,:)) %t-test %distractor absent
% [~,p] = ttest(sbj_meanrt_s2_da(1,:), sbj_meanrt_s2_da(2,:)) %t-test
% [~,p] = ttest(sbj_meanrt_ov_da(1,:), sbj_meanrt_ov_da(2,:)) %t-test
% [~,p] = ttest(sbj_meanrt_s1_dp(1,:), sbj_meanrt_s1_dp(2,:)) %t-test %distractor present
% [~,p] = ttest(sbj_meanrt_s2_dp(1,:), sbj_meanrt_s2_dp(2,:)) %t-test
% [~,p] = ttest(sbj_meanrt_ov_dp(1,:), sbj_meanrt_ov_dp(2,:)) %t-test

figure; subplot(231); boxplot(sbj_meanrt_s1_da'); hold on; plot(mean(sbj_meanrt_s1_da,2)','db'); %distractor absent
ylim([350, 650]); ylabel('RTs'); xticklabels({'high','low'}); title('Session 1 - Absent');
subplot(232); boxplot(sbj_meanrt_s2_da'); hold on; plot(mean(sbj_meanrt_s2_da,2)','db');
ylim([350, 650]); xticklabels({'high','low'}); title('Session 2- Absent');
subplot(233); boxplot(sbj_meanrt_ov_da'); hold on; plot(mean(sbj_meanrt_ov_da,2)','db');
ylim([350, 650]); xticklabels({'high','low'}); title('Overall- Absent');
subplot(234); boxplot(sbj_meanrt_s1_dp'); hold on; plot(mean(sbj_meanrt_s1_dp,2)','db'); %distractor present
ylim([350, 650]); ylabel('RTs'); xticklabels({'high','low'}); title('Session 1 - Present');
subplot(235); boxplot(sbj_meanrt_s2_dp'); hold on; plot(mean(sbj_meanrt_s2_dp,2)','db');
ylim([350, 650]); xticklabels({'high','low'}); title('Session 2 - Present');
subplot(236); boxplot(sbj_meanrt_ov_dp'); hold on; plot(mean(sbj_meanrt_ov_dp,2)','db');
ylim([350, 650]); xticklabels({'high','low'}); title('Overall - Present');


%% Target location X Distractor Presence - Accuracy analysis
% Mean accuracy
[means,stds,grp] = grpstats(data.correct,...
    {data.subject_nr, data.session_count, data.target_prob, data.distractor_pres},...
    {'mean','std','gname'});
means2 = reshape(means,2,2,2,sbj_nr);
stds2 = reshape(stds,2,2,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3),grp(:,4)),2,2,2,sbj_nr);


sbj_meanacc_s1_da = squeeze(means2(1,:,1,:)); %distractor absent
sbj_meanacc_s2_da = squeeze(means2(1,:,2,:));
sbj_meanacc_ov_da = squeeze(mean(means2(1,:,:,:),3));
sbj_meanacc_s1_dp = squeeze(means2(2,:,1,:)); %distractor present
sbj_meanacc_s2_dp = squeeze(means2(2,:,2,:));
sbj_meanacc_ov_dp = squeeze(mean(means2(2,:,:,:),3));

% [~,p] = ttest(sbj_meanacc_s1_da(1,:), sbj_meanacc_s1_da(2,:)) %t-test %distractor absent
% [~,p] = ttest(sbj_meanacc_s2_da(1,:), sbj_meanacc_s2_da(2,:)) %t-test
% [~,p] = ttest(sbj_meanacc_ov_da(1,:), sbj_meanacc_ov_da(2,:)) %t-test
% [~,p] = ttest(sbj_meanacc_s1_dp(1,:), sbj_meanacc_s1_dp(2,:)) %t-test %distractor present
% [~,p] = ttest(sbj_meanacc_s2_dp(1,:), sbj_meanacc_s2_dp(2,:)) %t-test
% [~,p] = ttest(sbj_meanacc_ov_dp(1,:), sbj_meanacc_ov_dp(2,:)) %t-test
figure; subplot(231); boxplot(sbj_meanacc_s1_da'); hold on; plot(mean(sbj_meanacc_s1_da,2)','db'); %distractor absent
ylim([.8, 1]); ylabel('Accuracy'); xticklabels({'high','low'}); title('Session 1 - Absent');
subplot(232); boxplot(sbj_meanacc_s2_da'); hold on; plot(mean(sbj_meanacc_s2_da,2)','db');
ylim([.8, 1]); xticklabels({'high','low'}); title('Session 2- Absent');
subplot(233); boxplot(sbj_meanacc_ov_da'); hold on; plot(mean(sbj_meanacc_ov_da,2)','db');
ylim([.8, 1]); xticklabels({'high','low'}); title('Overall- Absent');
subplot(234); boxplot(sbj_meanacc_s1_dp'); hold on; plot(mean(sbj_meanacc_s1_dp,2)','db'); %distractor present
ylim([.8, 1]); ylabel('Accuracy'); xticklabels({'high','low'}); title('Session 1 - Present');
subplot(235); boxplot(sbj_meanacc_s2_dp'); hold on; plot(mean(sbj_meanacc_s2_dp,2)','db');
ylim([.8, 1]); xticklabels({'high','low'}); title('Session 2 - Present');
subplot(236); boxplot(sbj_meanacc_ov_dp'); hold on; plot(mean(sbj_meanacc_ov_dp,2)','db');
ylim([.8, 1]); xticklabels({'high','low'}); title('Overall - Present');


%% Target location X Distractor Presence - IES analysis

sbj_meanies_s1_da = sbj_meanrt_s1_da ./ sbj_meanacc_s1_da;
sbj_meanies_s2_da = sbj_meanrt_s2_da ./ sbj_meanacc_s2_da;
sbj_meanies_ov_da = mean(cat(3,sbj_meanies_s1_da, sbj_meanies_s2_da),3);
sbj_meanies_s1_dp = sbj_meanrt_s1_dp ./ sbj_meanacc_s1_dp;
sbj_meanies_s2_dp = sbj_meanrt_s2_dp ./ sbj_meanacc_s2_dp;
sbj_meanies_ov_dp = mean(cat(3,sbj_meanies_s1_dp, sbj_meanies_s2_dp),3);

% [~,p] = ttest(sbj_meanies_s1_da(1,:), sbj_meanies_s1_da(2,:)) %t-test
% [~,p] = ttest(sbj_meanies_s2_da(1,:), sbj_meanies_s2_da(2,:)) %t-test
% [~,p] = ttest(sbj_meanies_ov_da(1,:), sbj_meanies_ov_da(2,:)) %t-test
% [~,p] = ttest(sbj_meanies_s1_dp(1,:), sbj_meanies_s1_dp(2,:)) %t-test
% [~,p] = ttest(sbj_meanies_s2_dp(1,:), sbj_meanies_s2_dp(2,:)) %t-test
% [~,p] = ttest(sbj_meanies_ov_dp(1,:), sbj_meanies_ov_dp(2,:)) %t-test

figure; subplot(231); boxplot(sbj_meanies_s1_da'); hold on; plot(mean(sbj_meanies_s1_da,2)','db');
ylim([400, 700]); ylabel('IES'); xticklabels({'high','low'}); title('Session 1 - Absent');
subplot(232); boxplot(sbj_meanies_s2_da'); hold on; plot(mean(sbj_meanies_s2_da,2)','db');
ylim([400, 700]); xticklabels({'high','low'}); title('Session 2 - Absent'); %xlabel('Target location');
subplot(233); boxplot(sbj_meanies_ov_da'); hold on; plot(mean(sbj_meanies_ov_da,2)','db');
ylim([400, 700]); xticklabels({'high','low'}); title('Overall - Absent');
subplot(234); boxplot(sbj_meanies_s1_dp'); hold on; plot(mean(sbj_meanies_s1_dp,2)','db');
ylim([400, 700]); ylabel('IES'); xticklabels({'high','low'}); title('Session 1 - Present');
subplot(235); boxplot(sbj_meanies_s2_dp'); hold on; plot(mean(sbj_meanies_s2_dp,2)','db');
ylim([400, 700]); xticklabels({'high','low'}); title('Session 2 - Present'); %xlabel('Target location');
subplot(236); boxplot(sbj_meanies_ov_dp'); hold on; plot(mean(sbj_meanies_ov_dp,2)','db');
ylim([400, 700]); xticklabels({'high','low'}); title('Overall - Present');

%% Frequency-tagged distractor location
% Mean reaction time
[means,stds,grp] = grpstats(datart.RT,...
    {datart.subject_nr, datart.session_count, datart.tagged_disloc},...
    {'mean','std','gname'});
means2 = reshape(means,3,2,sbj_nr);
stds2 = reshape(stds,3,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),3,2,sbj_nr);

sbj_meanrt_ov = mean(means2,2);
sbj_stdrt = mean(stds2,2);

group_meanrt = mean(sbj_meanrt_ov,3);
group_stdrt = mean(sbj_stdrt,3);
group_stert = group_stdrt/sqrt(sbj_nr);

sbj_meanrt_s1 = means2(:,1,:);
sbj_meanrt_s2 = means2(:,2,:);
sbj_meanrt_ov = mean(means2,2);
sbj_stdrt_s1 = stds2(:,1,:);
sbj_stdrt_s2 = stds2(:,2,:);

% [~,p] = ttest(sbj_meanrt_s1(2,:), sbj_meanrt_s1(3,:)) %t-test
% [~,p] = ttest(sbj_meanrt_s2(2,:), sbj_meanrt_s2(3,:)) %t-test
% [~,p] = ttest(sbj_meanrt_ov(2,:), sbj_meanrt_ov(3,:)) %t-test

figure; subplot(231); boxplot(squeeze(sbj_meanrt_s1)'); hold on; plot(mean(sbj_meanrt_s1,3)','db');
ylim([350, 650]); ylabel('RTs'); xticklabels({'absent','yes','no'}); title('Session 1');
subplot(232); boxplot(squeeze(sbj_meanrt_s2)'); hold on; plot(mean(sbj_meanrt_s2,3)','db');
ylim([350, 650]); xticklabels({'absent','yes','no'}); title('Session 2'); %xlabel('Distractor location');
subplot(233); boxplot(squeeze(sbj_meanrt_ov)'); hold on; plot(mean(sbj_meanrt_ov,3)','db');
ylim([350, 650]); xticklabels({'absent','yes','no'}); title('Overall');

% group_meanrt_s1 = mean(sbj_meanrt_s1,3);
% group_meanrt_s2 = mean(sbj_meanrt_s2,3);
% group_stdrt_s1 = mean(sbj_stdrt_s1,3);
% group_stdrt_s2 = mean(sbj_stdrt_s2,3);
% group_stert_s1 = group_stdrt_s1/sqrt(sbj_nr);
% group_stert_s2 = group_stdrt_s2/sqrt(sbj_nr);
% 
% figure; subplot(231); bar([group_meanrt_s1(1),group_meanrt_s1(3),group_meanrt_s1(2)]);
% ylim([350, 650]); ylabel('RTs'); xticklabels({'absent','no','yes'}); title('Session 1');
% hold on;
% errorbar([group_meanrt_s1(1),group_meanrt_s1(3),group_meanrt_s1(2)],[group_stert_s1(1),group_stert_s1(3),group_stert_s1(2)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(232); bar([group_meanrt_s2(1),group_meanrt_s2(3),group_meanrt_s2(2)]);
% ylim([350, 650]); xlabel('Frequency-tagged distractor location'); xticklabels({'absent','no','yes'}); title('Session 2'); %xlabel('Frequency-tagged distractor location');
% hold on;
% errorbar([group_meanrt_s2(1),group_meanrt_s2(3),group_meanrt_s2(2)],[group_stert_s2(1),group_stert_s2(3),group_stert_s2(2)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(233); bar([group_meanrt(1),group_meanrt(3),group_meanrt(2)]);
% ylim([350, 650]); xticklabels({'absent','no','yes'}); title('Overall');
% hold on;
% errorbar([group_meanrt(1),group_meanrt(3),group_meanrt(2)],[group_stert(1),group_stert(3),group_stert(2)],...
%     'LineStyle','none', 'Color', 'k');


%% Frequency-tagged target location
% Mean reaction time
[means,stds,grp] = grpstats(datart.RT,...
    {datart.subject_nr, datart.session_count, datart.tagged_tarloc},...
    {'mean','std','gname'});
means2 = reshape(means,2,2,sbj_nr);
stds2 = reshape(stds,2,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),2,2,sbj_nr);

sbj_meanrt_ov = mean(means2,2);
sbj_stdrt = mean(stds2,2);

group_meanrt = mean(sbj_meanrt_ov,3);
group_stdrt = mean(sbj_stdrt,3);
group_stert = group_stdrt/sqrt(sbj_nr);

sbj_meanrt_s1 = means2(:,1,:);
sbj_meanrt_s2 = means2(:,2,:);
sbj_meanrt_ov = mean(means2,2);
sbj_stdrt_s1 = stds2(:,1,:);
sbj_stdrt_s2 = stds2(:,2,:);

% [~,p] = ttest(sbj_meanrt_s1(1,:), sbj_meanrt_s1(2,:)) %t-test
% [~,p] = ttest(sbj_meanrt_s2(1,:), sbj_meanrt_s2(2,:)) %t-test
% [~,p] = ttest(sbj_meanrt_ov(1,:), sbj_meanrt_ov(2,:)) %t-test

subplot(234); boxplot(squeeze(sbj_meanrt_s1)');
ylim([350, 650]); ylabel('RTs'); xticklabels({'no','yes'}); title('Session 1');
subplot(235); boxplot(squeeze(sbj_meanrt_s2)');
ylim([350, 650]); xticklabels({'no','yes'}); title('Session 2'); %xlabel('Target location');
subplot(236); boxplot(squeeze(sbj_meanrt_ov)');
ylim([350, 650]); xticklabels({'no','yes'}); title('Overall');

% group_meanrt_s1 = mean(sbj_meanrt_s1,3);
% group_meanrt_s2 = mean(sbj_meanrt_s2,3);
% group_stdrt_s1 = mean(sbj_stdrt_s1,3);
% group_stdrt_s2 = mean(sbj_stdrt_s2,3);
% group_stert_s1 = group_stdrt_s1/sqrt(sbj_nr);
% group_stert_s2 = group_stdrt_s2/sqrt(sbj_nr);
% 
% subplot(234); bar([group_meanrt_s1(1),group_meanrt_s1(2)]);
% ylim([350, 650]); ylabel('RTs'); xticklabels({'no','yes'}); %title('Session 1');
% hold on;
% errorbar([group_meanrt_s1(1),group_meanrt_s1(2)],[group_stert_s1(1),group_stert_s1(2)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(235); bar([group_meanrt_s2(1),group_meanrt_s2(2)]);
% ylim([350, 650]); xlabel('Frequency-tagged target location'); xticklabels({'no','yes'}); %title('Session 2');
% hold on;
% errorbar([group_meanrt_s2(1),group_meanrt_s2(2)],[group_stert_s2(1),group_stert_s2(2)],...
%     'LineStyle','none', 'Color', 'k');
% subplot(236); bar([group_meanrt(1),group_meanrt(2)]);
% ylim([350, 650]); xticklabels({'no','yes'}); %title('Overall');
% hold on;
% errorbar([group_meanrt(1),group_meanrt(2)],[group_stert(1),group_stert(2)],...
%     'LineStyle','none', 'Color', 'k');


%% Absolute location - Distractor
% Mean reaction time
[means,stds,grp] = grpstats(datart.RT,...
    {datart.subject_nr, datart.session_count, datart.distractor_loc},...
    {'mean','std','gname'});
means2 = reshape(means,5,2,sbj_nr);
stds2 = reshape(stds,5,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),5,2,sbj_nr);

sbj_s1_rt = squeeze(means2(:,1,:));
sbj_s2_rt = squeeze(means2(:,2,:));
sbj_ov_rt = squeeze(mean(means2,2));

group_s1_rt = mean(sbj_s1_rt,2);
group_s2_rt = mean(sbj_s2_rt,2);
group_ov_rt = mean(sbj_ov_rt,2);

figure; subplot(231); plot(sbj_s1_rt,'.'); hold on; plot(group_s1_rt, 'k','LineWidth',1.5);
xlim([.5 5.5]); xticks(1:5); ylim([350, 700]); ylabel('RTs'); xticklabels({'ABS','tL','bL','tR','bR'}); title('Session 1');
subplot(232); plot(sbj_s2_rt,'.'); hold on; plot(group_s2_rt, 'k','LineWidth',1.5);
xlim([.5 5.5]); xticks(1:5); ylim([350, 700]); ylabel('RTs'); xticklabels({'ABS','tL','bL','tR','bR'}); title('Session 2');
subplot(233); plot(sbj_ov_rt,'.'); hold on; plot(group_ov_rt, 'k','LineWidth',1.5);
xlim([.5 5.5]); xticks(1:5); ylim([350, 700]); ylabel('RTs'); xticklabels({'ABS','tL','bL','tR','bR'}); title('Overall');

% Mean accuracy
[means,stds,grp] = grpstats(data.correct,...
    {data.subject_nr, data.session_count, data.distractor_loc},...
    {'mean','std','gname'});
means2 = reshape(means,5,2,sbj_nr);
stds2 = reshape(stds,5,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),5,2,sbj_nr);

sbj_s1_acc = squeeze(means2(:,1,:));
sbj_s2_acc = squeeze(means2(:,2,:));
sbj_ov_acc = squeeze(mean(means2,2));

group_s1_acc = mean(sbj_s1_acc,2);
group_s2_acc = mean(sbj_s2_acc,2);
group_ov_acc = mean(sbj_ov_acc,2);

subplot(234); plot(sbj_s1_acc,'.'); hold on;  plot(group_s1_acc, 'k','LineWidth',1.5);
xlim([.5 5.5]); xticks(1:5); ylim([.6, 1]); ylabel('Accuracy'); xticklabels({'ABS','tL','bL','tR','bR'}); title('Session 1');
subplot(235); plot(sbj_s2_acc,'.'); hold on; plot(group_s2_acc, 'k','LineWidth',1.5);
xlim([.5 5.5]); xticks(1:5); ylim([.6, 1]); ylabel('Accuracy'); xticklabels({'ABS','tL','bL','tR','bR'}); title('Session 2');
subplot(236); plot(sbj_ov_acc,'.'); hold on; plot(group_ov_acc, 'k','LineWidth',1.5);
xlim([.5 5.5]); xticks(1:5); ylim([.6, 1]); ylabel('Accuracy'); xticklabels({'ABS','tL','bL','tR','bR'}); title('Overall');


%% Absolute location - Target
% Mean reaction time
[means,stds,grp] = grpstats(datart.RT,...
    {datart.subject_nr, datart.session_count, datart.target_loc},...
    {'mean','std','gname'});
means2 = reshape(means,4,2,sbj_nr);
stds2 = reshape(stds,4,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),4,2,sbj_nr);

sbj_s1_rt = squeeze(means2(:,1,:));
sbj_s2_rt = squeeze(means2(:,2,:));
sbj_ov_rt = squeeze(mean(means2,2));

group_s1_rt = mean(sbj_s1_rt,2);
group_s2_rt = mean(sbj_s2_rt,2);
group_ov_rt = mean(sbj_ov_rt,2);

figure; subplot(231); plot(sbj_s1_rt,'.'); hold on; plot(group_s1_rt, 'k','LineWidth',1.5);
xlim([.5 4.5]); xticks(1:4); ylim([350, 700]); ylabel('RTs'); xticklabels({'tL','bL','tR','bR'}); title('Session 1');
subplot(232); plot(sbj_s2_rt,'.'); hold on; plot(group_s2_rt, 'k','LineWidth',1.5);
xlim([.5 4.5]); xticks(1:4); ylim([350, 700]); ylabel('RTs'); xticklabels({'tL','bL','tR','bR'}); title('Session 2');
subplot(233); plot(sbj_ov_rt,'.'); hold on; plot(group_ov_rt, 'k','LineWidth',1.5);
xlim([.5 4.5]); xticks(1:4); ylim([350, 700]); ylabel('RTs'); xticklabels({'tL','bL','tR','bR'}); title('Overall');

% Mean accuracy
[means,stds,grp] = grpstats(data.correct,...
    {data.subject_nr, data.session_count, data.target_loc},...
    {'mean','std','gname'});
means2 = reshape(means,4,2,sbj_nr);
stds2 = reshape(stds,4,2,sbj_nr);
grp2 = reshape(strcat(grp(:,1),grp(:,2),grp(:,3)),4,2,sbj_nr);

sbj_s1_acc = squeeze(means2(:,1,:));
sbj_s2_acc = squeeze(means2(:,2,:));
sbj_ov_acc = squeeze(mean(means2,2));

group_s1_acc = mean(sbj_s1_acc,2);
group_s2_acc = mean(sbj_s2_acc,2);
group_ov_acc = mean(sbj_ov_acc,2);

subplot(234); plot(sbj_s1_acc,'.'); hold on;  plot(group_s1_acc, 'k','LineWidth',1.5);
xlim([.5 4.5]); xticks(1:4); ylim([.7, 1]); ylabel('Accuracy'); xticklabels({'tL','bL','tR','bR'}); title('Session 1');
subplot(235); plot(sbj_s2_acc,'.'); hold on; plot(group_s2_acc, 'k','LineWidth',1.5);
xlim([.5 4.5]); xticks(1:4); ylim([.7, 1]); ylabel('Accuracy'); xticklabels({'tL','bL','tR','bR'}); title('Session 2');
subplot(236); plot(sbj_ov_acc,'.'); hold on; plot(group_ov_acc, 'k','LineWidth',1.5);
xlim([.5 4.5]); xticks(1:4); ylim([.7, 1]); ylabel('Accuracy'); xticklabels({'tL','bL','tR','bR'}); title('Overall');

