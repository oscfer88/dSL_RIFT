% Cross-correlation analysis

ft_defaults

clearvars

%% Select data files and event types
%sbj1_session1
cd('Q:/MEG_data/20190724_b64d/190724/');
dataset{1,1} = 'B64D_S1P1.fif';
dataset{1,2} = 'B64D_S1P2.fif';
%sbj1_session2
cd('Q:/MEG_data/20190801_b64d/190801/');
% dataset{1,1} = 'B64D_S2P1.fif';
dataset{1,1} = 'B64D_S2P2.fif'; %dataset{1,2}
%sbj2_session1
cd('Q:/MEG_data/20190726_b32e/190726/');
dataset{1,1} = 'B32E_S1P1.fif';
dataset{1,2} = 'B32E_S1P2.fif';
%sbj2_session2
cd('Q:/MEG_data/20190801_b32e/190801/');
dataset{1,1} = 'B32E_S2P1.fif';
dataset{1,2} = 'B32E_S2P2.fif';
%sbj3_session1
cd('Q:/MEG_data/20190726_b581/190726/');
dataset{1,1} = 'B581_S1P1.fif';
dataset{1,2} = 'B581_S1P2.fif';
%sbj3_session2
cd('Q:/MEG_data/20190801_b581/190801/');
dataset{1,1} = 'B581_S2P1.fif';
dataset{1,2} = 'B581_S2P2.fif';
%sbj4_session1
cd('Q:/MEG_data/20190730_b32d/190730/');
dataset{1,1} = 'B32D_S1P1.fif';
dataset{1,2} = 'B32D_S1P2.fif';
%sbj4_session2
cd('Q:/MEG_data/20190805_b32d/190805/');
dataset{1,1} = 'B32D_S2P1.fif';
dataset{1,2} = 'B32D_S2P2.fif';
%sbj5_session1
cd('Q:/MEG_data/20190830_b3ed/190830/');
dataset{1,1} = 'B3ED_S1P1.fif';
dataset{1,2} = 'B3ED_S1P2.fif';
%sbj5_session2
cd('Q:/MEG_data/20190903_b3ed/190903/');
dataset{1,1} = 'B3ED_S2P1.fif';
dataset{1,2} = 'B3ED_S2P2.fif';
%sbj6_session1
cd('Q:/MEG_data/20190830_b4b8/190830/');
dataset{1,1} = 'B4B8_S1P1.fif';
dataset{1,2} = 'B4B8_S1P2.fif';
%sbj6_session2
cd('Q:/MEG_data/20190903_b4b8/190903/');
dataset{1,1} = 'B4B8_S2P1.fif';
dataset{1,2} = 'B4B8_S2P2.fif';
%sbj7_session1
cd('Q:/MEG_data/20190903_b3f3/190903/');
dataset{1,1} = 'B3F3_S1P1.fif';
dataset{1,2} = 'B3F3_S1P2.fif';
%sbj7_session2
cd('Q:/MEG_data/20190911_b3f3/190911/');
dataset{1,1} = 'B3F3_S2P1.fif';
dataset{1,2} = 'B3F3_S2P2.fif';
%sbj8_session1
cd('Q:/MEG_data/20190909_b3f5/190909/');
dataset{1,1} = 'B3F5_S1P1.fif';
dataset{1,2} = 'B3F5_S1P2.fif';
%sbj8_session2
cd('Q:/MEG_data/20190913_b3f5/190913/');
dataset{1,1} = 'B3F5_S2P1.fif';
dataset{1,2} = 'B3F5_S2P1-1.fif';
dataset{1,3} = 'B3F5_S2P2.fif';
%sbj9_session1
cd('Q:/MEG_data/20190916_b3f2/190916/');
dataset{1,1} = 'B3F2_S1P1.fif';
dataset{1,2} = 'B3F2_S1P2.fif';
%sbj9_session2
cd('Q:/MEG_data/20190920_b3f2/190920/');
dataset{1,1} = 'B3F2_S2P1.fif';
dataset{1,2} = 'B3F2_S2P1.fif';
%sbj10_session1
cd('Q:/MEG_data/20190916_b57f/190916/');
dataset{1,1} = 'B57F_S1P1.fif';
dataset{1,2} = 'B57F_S1P1.fif';
%sbj10_session2
cd('Q:/MEG_data/20190920_b57f/190920/');
dataset{1,1} = 'B57F_S2P1.fif';
dataset{1,2} = 'B57F_S2P1.fif';
%sbj11_session1
cd('Q:/MEG_data/20191014_b3f3/191014/');
dataset{1,1} = 'B3F3_S1P1.fif';
dataset{1,2} = 'B3F3_S1P2.fif';
%sbj11_session2
cd('Q:/MEG_data/20191015_b3f3/191015/');
dataset{1,1} = 'B3F3_S2P1.fif';
dataset{1,2} = 'B3F3_S2P2.fif';
%sbj12_session1
cd('Q:/MEG_data/20191014_b4c2/191014/');
dataset{1,1} = 'B4C2_S1P1.fif';
dataset{1,2} = 'B4C2_S1P2.fif';
%sbj12_session2
cd('Q:/MEG_data/20191015_b4c2/191015/');
dataset{1,1} = 'B4C2_S2P1.fif';
dataset{1,2} = 'B4C2_S2P2.fif';
%sbj13_session1
cd('Q:/MEG_data/20191107_b3a4/191107/');
dataset{1,1} = 'B3A4_S1P1.fif';
dataset{1,2} = 'B3A4_S1P2.fif';
%sbj13_session2
cd('Q:/MEG_data/20191113_b3a4/191113/');
dataset{1,1} = 'B3A4_S2P1.fif';
dataset{1,2} = 'B3A4_S2P2.fif';
%sbj14_session1
cd('Q:/MEG_data/20191111_b463/191111/');
dataset{1,1} = 'B463_S1P1.fif';
dataset{1,2} = 'B463_S1P2.fif';
%sbj14_session2
cd('Q:/MEG_data/20191113_b463/191113/');
dataset{1,1} = 'B463_S2P1.fif';
dataset{1,2} = 'B463_S2P2.fif';
%sbj15_session1
cd('Q:/MEG_data/20191111_b5f2/191111/');
dataset{1,1} = 'B5F2_S1P1.fif';
dataset{1,2} = 'B5F2_S1P2.fif';
%sbj15_session2
cd('Q:/MEG_data/20191115_b5f2/191115/');
dataset{1,1} = 'B5F2_S2P1.fif';
dataset{1,2} = 'B5F2_S2P2.fif';
%sbj16_session1
cd('Q:/MEG_data/20191115_b45f/191115/');
dataset{1,1} = 'B45F_S1P1.fif';
dataset{1,2} = 'B45F_S1P2.fif';
%sbj16_session2
cd('Q:/MEG_data/20191120_b45f/191120/');
dataset{1,1} = 'B45F_S2P1.fif';
dataset{1,2} = 'B45F_S2P2.fif';
%sbj17_session1
cd('Q:/MEG_data/20191119_b521/191119/');
dataset{1,1} = 'B521_S1P1.fif';
dataset{1,2} = 'B521_S1P2.fif';
%sbj17_session2
cd('Q:/MEG_data/20191120_b521/191120/');
dataset{1,1} = 'B521_S2P1.fif';
dataset{1,2} = 'B521_S2P2.fif';
%sbj18_session1
cd('Q:/MEG_data/20191119_b5e7/191119/');
dataset{1,1} = 'B5E7_S1P1.fif';
dataset{1,2} = 'B5E7_S1P2.fif';
%sbj18_session2
cd('Q:/MEG_data/20191120_b5e7/191120/');
dataset{1,1} = 'B5E7_S2P1.fif';
dataset{1,2} = 'B5E7_S2P2.fif';
%sbj19_session1
cd('Q:/MEG_data/20191125_b394/191125/');
dataset{1,1} = 'B394_S1P1.fif';
dataset{1,2} = 'B394_S1P2.fif';
%sbj19_session2
cd('Q:/MEG_data/20191127_b394/191127/');
dataset{1,1} = 'B394_S2P1.fif';
dataset{1,2} = 'B394_S2P2.fif';
%sbj20_session1
cd('Q:/MEG_data/20191125_b45f/191125/');
dataset{1,1} = 'B45F_S1P1.fif';
dataset{1,2} = 'B45F_S1P2.fif';
%sbj20_session2
cd('Q:/MEG_data/20191128_b45f/191128/');
dataset{1,1} = 'B45F_S2P1.fif';
dataset{1,2} = 'B45F_S2P2.fif';

%Get fequency tagging start trigger
trigger_eventtype = [2,10,18,26,34,42,50,58,66,74,82,90,98,106,114,122,130,138,146,154,162,170,178,186,194,202,210,218,226,234,242,250];

%% Reading and segmenting the data creating trials
for i = 1:size(dataset,2)

    cfg                         = [];
    cfg.dataset                 = char(dataset{1,i});
    cfg.trialfun                = 'ft_trialfun_general';  
    cfg.trialdef.eventtype      = 'STI101';  
    cfg.trialdef.eventvalue     = trigger_eventtype;    
    cfg.trialdef.prestim        = 1.5; 
    cfg.trialdef.poststim       = 2.5; 
    cfg = ft_definetrial(cfg);
    
%     cfg.lpfilter                = 'yes';
%     cfg.lpfreq                  = 100;
%     cfg.padding                 = 5;     
%     cfg.padtype                 = 'data';
    cfg.channel                 = {'MEGGRAD' 'MISC004' 'MISC005'};
    cfg.detrend                 = 'yes';    
    data(i,1) = ft_preprocessing(cfg); 

end

% Merge the datafiles
if size(data,1) == 1
    data_planar = ft_appenddata([], data(1,1));
elseif size(data,1) == 2
    data_planar = ft_appenddata([], data(1,1),data(2,1));
elseif size(data,1) == 3
    data_planar = ft_appenddata([], data(1,1),data(2,1),data(3,1));
end
data_planar = rmfield(data_planar,'sampleinfo'); %remove the 'sampleinfo' field of the newly created datafile

%% Artifact rejections in planar gradiometers
load trl_keep
load good_trials %correct response trials
load chan_keep
load eye_keep2
nolear_trials = 101:720; %remove learning (the first 100 trials)
cfg                     = [];
cfg.trials              = intersect(trl_keep, intersect(eye_keep, intersect(good_trials, nolear_trials)));
cfg.channel             = [chan_keep(:);  {'MISC004'}; {'MISC005'}];
planar_rjv = ft_selectdata(cfg,data_planar);

% Artifact suppression by ICA
cfg = [];
cfg.channel = chan_keep(:);
planar_rjv2 = ft_selectdata(cfg,planar_rjv);

load comp3
cfg = [];
cfg.component = [11 15 23];       %sbj1_s1
% cfg.component = [11 17 20];       %sbj1_s2
% cfg.component = [2 11 15];        %sbj2_s1
% cfg.component = [1 12 18];        %sbj2_s2
% cfg.component = [5 6 9];          %sbj3_s1
% cfg.component = [2 10 11];        %sbj3_s2
% cfg.component = [24 28];          %sbj4_s1
% cfg.component = [22 25];          %sbj4_s2
% cfg.component = [1 17 22];        %sbj5_s1
% cfg.component = [1 10];           %sbj5_s2
% cfg.component = [1 4 8];          %sbj6_s1
% cfg.component = [1 9];            %sbj6_s2
% cfg.component = [1 5 13];         %sbj7_s1
% cfg.component = [2 4 9];          %sbj7_s2
% cfg.component = [1 2];            %sbj8_s1
% cfg.component = [1 2];            %sbj8_s2
% cfg.component = [2 13];           %sbj9_s1
% cfg.component = [4 9];            %sbj9_s2
% cfg.component = [1 4];            %sbj10_s1
% cfg.component = [2 4];            %sbj10_s2
% cfg.component = [24];             %sbj11_s1
% cfg.component = [20 21];          %sbj11_s2
% cfg.component = [1 2 10];         %sbj12_s1
% cfg.component = [1 2 7];          %sbj12_s2
% cfg.component = [7];              %sbj13_s1
% cfg.component = [8];              %sbj13_s2
% cfg.component = [5 9 19];         %sbj14_s1
% cfg.component = [2 19];           %sbj14_s2
% cfg.component = [1 2 8];          %sbj15_s1
% cfg.component = [1 4];            %sbj15_s2
% cfg.component = [9];              %sbj16_s1
% cfg.component = [2 11];           %sbj16_s2
% cfg.component = [1 3 9];          %sbj17_s1
% cfg.component = [3 5];            %sbj17_s2
% cfg.component = [1 4];            %sbj18_s1
% cfg.component = [1 4];            %sbj18_s2
% cfg.component = [7 38];           %sbj19_s1
% cfg.component = [14 42 46];       %sbj19_s2
% cfg.component = [3 17];           %sbj20_s1
% cfg.component = [1 22];           %sbj20_s2
planar_clear  = ft_rejectcomponent(cfg, planar_comp, planar_rjv2);


%% Cross-correlation
%prepare the data
timeStart = -0.2;
timeEnd = 1.7;

cfg                     = [];
cfg.latency             = [timeStart timeEnd];
datameg_xcor = ft_selectdata(cfg,planar_clear);
% save('datameg_xcor', 'datameg_xcor');

cfg                     = [];
cfg.channel             = 'MISC005'; %for sbj2_s1 & sbj3_s1: 'MISC004'
cfg.latency             = [timeStart timeEnd];
dataphotoL_xcor = ft_selectdata(cfg,planar_rjv);  %for sbj3_s1: missing 'MISC004'
% save('dataphotoL_xcor', 'dataphotoL_xcor');

cfg                     = [];
cfg.channel             = 'MISC004'; %for sbj2_s1 & sbj3_s1: 'MISC005'
cfg.latency             = [timeStart timeEnd];
dataphotoR_xcor = ft_selectdata(cfg,planar_rjv);
% save('dataphotoR_xcor', 'dataphotoR_xcor');

datameg_xcor = cat(3,datameg_xcor.trial{:});
dataphotoL_xcor = cat(3,dataphotoL_xcor.trial{:});
dataphotoR_xcor = cat(3,dataphotoR_xcor.trial{:});

%run cross-correlation function
% load datameg_xcor
% load dataphotoL_xcor
% load dataphotoR_xcor
[xcorL, xcorR, pDelays] = cross_correlation(datameg_xcor, dataphotoL_xcor, dataphotoR_xcor, timeStart, timeEnd);
% save('xcorL3','xcorL')
% save('xcorR3','xcorR')


% %%%
% % Load psychtoolbox/matlab data and control the frequency tagging signal
% signalL = tag_list(1,:);
% signalL2 = cell(1,360);
% for i=1:numel(signalL)
%     signalL2{i} = reshape(signalL{i}, 1, []);
% end
% signalL3 = nan(1,864,360);
% for i=1:numel(signalL2)
%     signalL3(:,:,i) = signalL2{i};
% end
% signalL4 = signalL3(:,1:(864/1.8)*1.5,:);
% 
% dataphotoL_xcor_resample = nan(1,720,360);
% for i=1:360
%     dataphotoL_xcor_resample(:,:,i) = resample(dataphotoL_xcor(:,:,i),720,1500); %part1
% %     dataphotoL_xcor_resample(:,:,i) = resample(dataphotoL_xcor(:,:,360+i),720,1500); %part2
% end
% 
% signalR = tag_list(2,:);
% signalR2 = cell(1,360);
% for i=1:numel(signalR)
%     signalR2{i} = reshape(signalR{i}, 1, []);
% end
% signalR3 = nan(1,864,360);
% for i=1:numel(signalR2)
%     signalR3(:,:,i) = signalR2{i};
% end
% signalR4 = signalR3(:,1:(864/1.8)*1.5,:);
% 
% dataphotoR_xcor_resample = nan(1,720,360);
% for i=1:360
%     dataphotoR_xcor_resample(:,:,i) = resample(dataphotoR_xcor(:,:,i),720,1500); %part1
% %     dataphotoR_xcor_resample(:,:,i) = resample(dataphotoR_xcor(:,:,360+i),720,1500); %part2
% end
% 
% [xcorL, xcorR, pDelays] = cross_correlation(dataphotoL_xcor_resample, signalL4, signalR4, 360);
% figure;
% subplot(221); plot(pDelays, xcorL); ylim([-1 1]); xlabel('delays'); ylabel('correlation');
% subplot(223); plot(pDelays, abs(hilbert(xcorL'))); ylim([0 1]); xlabel('delays'); ylabel('correlation');
% subplot(222); plot(pDelays, xcorR); ylim([-1 1]); xlabel('delays'); ylabel('correlation');
% subplot(224); plot(pDelays, abs(hilbert(xcorR'))); ylim([0 1]); xlabel('delays'); ylabel('correlation');
% 
% [xcorL, xcorR, pDelays] = cross_correlation(dataphotoR_xcor_resample, signalL4, signalR4, 360);
% figure;
% subplot(221); plot(pDelays, xcorL); ylim([-1 1]); xlabel('delays'); ylabel('correlation');
% subplot(223); plot(pDelays, abs(hilbert(xcorL'))); ylim([0 1]); xlabel('delays'); ylabel('correlation');
% subplot(222); plot(pDelays, xcorR); ylim([-1 1]); xlabel('delays'); ylabel('correlation');
% subplot(224); plot(pDelays, abs(hilbert(xcorR'))); ylim([0 1]); xlabel('delays'); ylabel('correlation');
% %%%


%% Plots

% compute hilbert transform
% load xcorL3
% load xcorR3
pDelays = -200:1:200; %
xcorLhilb = abs(hilbert(xcorL'));
xcorRhilb = abs(hilbert(xcorR'));

%plot cross-correlation 1
figure;
subplot(221); plot(pDelays, xcorL); ylim([-.20 .20]); xlabel('delays'); ylabel('correlation');
subplot(223); plot(pDelays, abs(hilbert(xcorL'))); ylim([0 .20]); xlabel('delays'); ylabel('correlation');
subplot(222); plot(pDelays, xcorR); ylim([-.20 .20]); xlabel('delays'); ylabel('correlation');
subplot(224); plot(pDelays, abs(hilbert(xcorR'))); ylim([0 .20]); xlabel('delays'); ylabel('correlation');


% %plot cross-correlation 2
% figure;
% j = 0;
% for i=1:203
%     if i == 109
%         figure
%         j = 108;
%     end
%     idx = i-j;
%     subplot(9,12,idx);
%     plot(pDelays, abs(hilbert(xcorL(idx,:)')));
%     hold on;
%     plot(pDelays, abs(hilbert(xcorR(idx,:)')));
%     ylim([0 .2]); xlabel('delays'); ylabel('correlation');
% end
% legend({'left' 'right'});


% Multiplot cross-correlation
load tfr_all
tfr_all.freq = 1;
tfr_all.freq = pDelays;
tfr_all.dimord = 'chan_freq';
for i=1:2
    if i == 1
        tfr_all.powspctrm = xcorLhilb';
    elseif i == 2
        tfr_all.powspctrm = xcorRhilb';
    end
%     tfr_all.powspctrm = pcorLhilb' - pcorRhilb';
    %plot all planar sensors
    cfg                 = [];
    cfg.ylim            = [0 0.15];
%     cfg.ylim            = [-0.15 0.15];
    cfg.layout          = 'neuromag306planar.lay';
    figure;ft_multiplotER(cfg, tfr_all);
%     %combine sensors
%     cfg             = [];
%     cfg.method      = 'sum';
%     tfr_all_c = ft_combineplanar(cfg,tfr_all);
%     %topoplot
%     cfg                 = [];
%     cfg.ylim            = [0 .2] ;	        
%     cfg.layout          = 'neuromag306cmb.lay';
%     figure;ft_multiplotER(cfg,tfr_all_c);
end


% Topoplot cross-correlation
load tfr_all
tfr_all.freq = 1;
tfr_all.time = 1;
tfr_all.dimord = 'chan_freq';
for i=1:2
    if i == 1
        xcor = xcorLhilb';
    elseif i == 2
        xcor = xcorRhilb';
    end
    
    %get values
%     xcor_value = xcor(:,141);
    xcor_value = max(xcor,[],2);

    %replace values
    tfr_all.powspctrm = xcor_value;
    
%     %smoothing1
%     tfr_all.powspctrm = smooth(xcor_value);
%     %smoothing2
%     windowSize = 5; 
%     b = (1/windowSize)*ones(1,windowSize);
%     a = 1;
%     tfr_all.powspctrm = filter(b,a,xcor_value);

    %combine sensors
    cfg             = [];
    cfg.method      = 'sum';
    tfr_all_c = ft_combineplanar(cfg,tfr_all);
    
    %topoplot
    figure;
    cfg                 = [];
    cfg.zlim            = [0 .20] ;	        
    cfg.layout          = 'neuromag306cmb.lay';
    plot_mat = ft_topoplotER(cfg,tfr_all_c);
    colormap jet
end


% Topoplot left > right
load tfr_all
tfr_all.freq = 1;
tfr_all.time = 1;
tfr_all.dimord = 'chan_freq';
%get values
% xcorLhilb_value = xcorLhilb(141,:);
xcorLhilb_value = max(xcorLhilb,[],1);
% xcorRhilb_value = xcorRhilb(141,:);
xcorRhilb_value = max(xcorRhilb,[],1);
xcorLmR = (xcorLhilb_value-xcorRhilb_value);%./(xcorLhilb_value+xcorRhilb_value);
%replace values
tfr_all.powspctrm = xcorLmR';
%combine sensors
cfg             = [];
cfg.method      = 'sum';
tfr_all_c = ft_combineplanar(cfg,tfr_all);
%topoplot
cfg = [];
cfg.zlim            = [-.1 .1] ;
% cfg.zlim            = [-.1 .1] ;
cfg.layout = 'neuromag306cmb.lay';
ft_topoplotER(cfg,tfr_all_c)
colormap jet


% Topoplot as a function of time
load tfr_all
tfr_all = rmfield(tfr_all,'freq');
tfr_all.time = pDelays/1000;
tfr_all.dimord = 'chan_time';
for i=1:2
    %values
    if i == 1
        tfr_all.powspctrm = xcorLhilb';
    elseif i == 2
        tfr_all.powspctrm = xcorRhilb';
    end
    %combine sensors
    cfg             = [];
    cfg.method      = 'sum';
    tfr_all_c = ft_combineplanar(cfg,tfr_all);
    %topoplot
    timepoint = [-.125:.025:.025];
    figure(i);
    for ii=1:length(timepoint)-1
        cfg                 = []; %plot topography of alpha (8-12 HZ) as a function of time
        cfg.xlim            = [timepoint(ii) timepoint(ii+1)];
        cfg.zlim            = [0 .15];
        cfg.layout          = 'neuromag306cmb.lay';
        subplot(2,3,ii);ft_topoplotTFR(cfg, tfr_all_c);
        colormap jet
    end
end


% Peak/time analysis
idx_pDelays_m100_0 = find(pDelays == -100):find(pDelays == 0);
pDelays_m100_0 = pDelays(idx_pDelays_m100_0);
xcorLhilbX = xcorLhilb(idx_pDelays_m100_0,:);
xcorRhilbX = xcorRhilb(idx_pDelays_m100_0,:);

% figure;
% subplot(121);
% plot(pDelays_m100_0,xcorLhilbX); ylim([0 .2]); xlabel('delays'); ylabel('correlation');
% subplot(122);
% plot(pDelays_m100_0,xcorRhilbX); ylim([0 .2]); xlabel('delays'); ylabel('correlation');

for i=1:size(xcorLhilbX,2)
    if max(xcorLhilbX(:,i)) < .02
        xcorLhilbX(:,i) = 0;
    end
    if max(xcorRhilbX(:,i)) < .02
        xcorRhilbX(:,i) = 0;
    end
end

% figure;
% subplot(121);
% plot(pDelays_m100_0,xcorLhilbX); ylim([0 .2]); xlabel('delays'); ylabel('correlation');
% subplot(122);
% plot(pDelays_m100_0,xcorRhilbX); ylim([0 .2]); xlabel('delays'); ylabel('correlation');

[xcorLhilbX_M,xcorLhilbX_I] = max(xcorLhilbX);
[xcorRhilbX_M,xcorRhilbX_I] = max(xcorRhilbX);

xcorLhilbX_D = pDelays_m100_0(xcorLhilbX_I);
xcorRhilbX_D = pDelays_m100_0(xcorRhilbX_I);

load tfr_all
tfr_all.freq = 1;
tfr_all.time = 1;
tfr_all.dimord = 'chan_freq';
for i=1:2
    if i == 1
        tfr_all.powspctrm = abs(xcorLhilbX_D)';
    elseif i == 2
        tfr_all.powspctrm = abs(xcorRhilbX_D)';
    end
    %combine sensors
    cfg             = [];
    cfg.method      = 'sum';
    tfr_all_c = ft_combineplanar(cfg,tfr_all);
    tfr_all_c.powspctrm = tfr_all_c.powspctrm/2;
    %topoplot
    figure;
    cfg                 = [];
    cfg.zlim            = [50 100] ;	        
    cfg.layout          = 'neuromag306cmb.lay';
    ft_topoplotER(cfg,tfr_all_c);
    colormap(flipud(jet));colorbar
end

% %% Partial correlation
% [pcorL, pcorR, ~] = partial_correlation(datameg_xcor, dataphotoL_xcor, dataphotoR_xcor, nTrials);
% % save('pcorL','pcorL')
% % save('pcorR','pcorR')
% 
% % compute hilbert transform
% pcorLhilb = abs(hilbert(pcorL'));
% pcorRhilb = abs(hilbert(pcorR'));
% pDelays = -200:1:200;
% 
% %plot partial correlation 1
% % load pcorL
% % load pcorR
% figure;
% subplot(221); plot(pDelays, pcorL); ylim([-.2 .2]); xlabel('delays'); ylabel('correlation');
% subplot(223); plot(pDelays, pcorLhilb); ylim([0 .2]); xlabel('delays'); ylabel('correlation');
% subplot(222); plot(pDelays, pcorR); ylim([-.2 .2]); xlabel('delays'); ylabel('correlation');
% subplot(224); plot(pDelays, pcorRhilb); ylim([0 .2]); xlabel('delays'); ylabel('correlation');
% 
% 
% % Multiplot partial correlation
% load tfr_all
% tfr_all.freq = 1;
% tfr_all.freq = pDelays;
% tfr_all.dimord = 'chan_freq';
% for i=1:2
%     if i == 1
%         tfr_all.powspctrm = pcorLhilb';
%     elseif i == 2
%         tfr_all.powspctrm = pcorRhilb';
%     end
% %     tfr_all.powspctrm = pcorLhilb' - pcorRhilb';
%     cfg                 = []; %plot all planar sensors
%     cfg.ylim            = [0 0.15];
% %     cfg.ylim            = [-0.15 0.15];
%     cfg.layout          = 'neuromag306planar.lay';
%     figure;ft_multiplotER(cfg, tfr_all);
% %     %combine sensors
% %     cfg             = [];
% %     cfg.method      = 'sum';
% %     tfr_all_c = ft_combineplanar(cfg,tfr_all);
% %     %topoplot
% %     cfg                 = [];
% %     cfg.ylim            = [0 .2] ;	        
% %     cfg.layout          = 'neuromag306cmb.lay';
% %     figure;ft_multiplotER(cfg,tfr_all_c);
% end\
% 
% 
% % Topoplot partial correlation
% load tfr_all
% tfr_all.freq = 1;
% tfr_all.time = 1;
% tfr_all.dimord = 'chan_freq';
% for i=1:2
%     if i == 1
%         pcor = pcorLhilb';
%     elseif i == 2
%         pcor = pcorRhilb';
%     end
%     %get values
% %     pcor_value = pcor(:,141);
%     pcor_value = max(pcor,[],2);
%     %replace values
%     tfr_all.powspctrm = pcor_value;
%     %combine sensors
%     cfg             = [];
%     cfg.method      = 'sum';
%     tfr_all_c = ft_combineplanar(cfg,tfr_all);
%     %topoplot
%     figure;
%     cfg                 = [];
%     cfg.zlim            = [0 .2] ;	        
%     cfg.layout          = 'neuromag306cmb.lay';
%     ft_topoplotER(cfg,tfr_all_c);
%     colormap jet
% end
% 
% 
% % Topoplot as a function of time
% load tfr_all
% tfr_all = rmfield(tfr_all,'freq');
% tfr_all.time = pDelays/1000;
% tfr_all.dimord = 'chan_time';
% for i=1:2
%     %values
%     if i == 1
%         tfr_all.powspctrm = pcorLhilb';
%     elseif i == 2
%         tfr_all.powspctrm = pcorRhilb';
%     end
%     %combine sensors
%     cfg             = [];
%     cfg.method      = 'sum';
%     tfr_all_c = ft_combineplanar(cfg,tfr_all);
%     %topoplot
%     timepoint = [-.125:.025:.025];
%     figure(i);
%     for ii=1:length(timepoint)-1
%         cfg                 = []; %plot topography of alpha (8-12 HZ) as a function of time
%         cfg.xlim            = [timepoint(ii) timepoint(ii+1)];
%         cfg.zlim            = [0 .15];
%         cfg.layout          = 'neuromag306cmb.lay';
%         subplot(2,3,ii);ft_topoplotTFR(cfg, tfr_all_c);
%         colormap jet
%     end
% end
% 
% 
% % Peak/time analysis
% pcorLhilbX = pcorLhilb(find(pDelays == -100):find(pDelays == 0),:);
% pcorRhilbX = pcorRhilb(find(pDelays == -100):find(pDelays == 0),:);
% figure;subplot(121);plot(pDelays(find(pDelays == -100):find(pDelays == 0)),pcorLhilbX); ylim([0 .2]); xlabel('delays'); ylabel('correlation');
% subplot(122);plot(pDelays(find(pDelays == -100):find(pDelays == 0)),pcorRhilbX); ylim([0 .2]); xlabel('delays'); ylabel('correlation');
% 
% for i=1:size(pcorLhilbX,2)
%     if max(pcorLhilbX(:,i)) < .02
%         pcorLhilbX(:,i) = 0;
%     end
%     if max(pcorRhilbX(:,i)) < .02
%         pcorRhilbX(:,i) = 0;
%     end
% end
% 
% pDelays_m100_0 = pDelays(find(pDelays == -100):find(pDelays == 0));
% 
% figure;
% subplot(121);
% plot(pDelays_m100_0,pcorLhilbX); ylim([0 .2]); xlabel('delays'); ylabel('correlation');
% subplot(122);
% plot(pDelays_m100_0,pcorRhilbX); ylim([0 .2]); xlabel('delays'); ylabel('correlation');
% 
% [pcorLhilbX_M,pcorLhilbX_I] = max(pcorLhilbX);
% [pcorRhilbX_M,pcorRhilbX_I] = max(pcorRhilbX);
% 
% pcorLhilbX_D = pDelays_m100_0(pcorLhilbX_I);
% pcorRhilbX_D = pDelays_m100_0(pcorRhilbX_I);
% 
% load tfr_all
% tfr_all.freq = 1;
% tfr_all.time = 1;
% tfr_all.dimord = 'chan_freq';
% for i=1:2
%     if i == 1
%         tfr_all.powspctrm = abs(pcorLhilbX_D)';
%     elseif i == 2
%         tfr_all.powspctrm = abs(pcorRhilbX_D)';
%     end
%     %combine sensors
%     cfg             = [];
%     cfg.method      = 'sum';
%     tfr_all_c = ft_combineplanar(cfg,tfr_all);
%     tfr_all_c.powspctrm = tfr_all_c.powspctrm/2;
%     %topoplot
%     figure;
%     cfg                 = [];
%     cfg.zlim            = [50 100] ;	        
%     cfg.layout          = 'neuromag306cmb.lay';
%     ft_topoplotER(cfg,tfr_all_c);
%     colormap(flipud(jet));colorbar
% end
% 
% 
% % Topoplot left > right
% load tfr_all
% tfr_all.freq = 1;
% tfr_all.time = 1;
% tfr_all.dimord = 'chan_freq';
% %get values
% pcorLhilb_value = pcorLhilb(141,:);
% % pcorLhilb_value = max(pcorLhilb,[],1);
% pcorRhilb_value = pcorRhilb(141,:);
% % pcorRhilb_value = max(pcorRhilb,[],1);
% pcorLmR = (pcorLhilb_value-pcorRhilb_value);%./(pcorLhilb_value+pcorRhilb_value);
% %replace values
% tfr_all.powspctrm = pcorLmR';
% %combine sensors
% cfg             = [];
% cfg.method      = 'sum';
% tfr_all_c = ft_combineplanar(cfg,tfr_all);
% %topoplot
% cfg = [];
% cfg.zlim            = [-1 1] ;
% % cfg.zlim            = [-.1 .1] ;
% cfg.layout = 'neuromag306cmb.lay';
% ft_topoplotER(cfg,tfr_all_c)
% colormap jet
% 
% 
% %% Analysis on half dataset at a time (cross-correlation only)
% 
% %divede data
% datameg_xcor_fsthalf = datameg_xcor(:,1:end/2,:);
% dataphotoL_xcor_fsthalf = dataphotoL_xcor(:,1:end/2,:);
% dataphotoR_xcor_fsthalf = dataphotoR_xcor(:,1:end/2,:);
% 
% datameg_xcor_sndhalf = datameg_xcor(:,1+end/2:end,:);
% dataphotoL_xcor_sndhalf = dataphotoL_xcor(:,1+end/2:end,:);
% dataphotoR_xcor_sndhalf = dataphotoR_xcor(:,1+end/2:end,:);
% 
% %run cross-correlation function
% [xcorL_fsthalf, xcorR_fsthalf, pDelays_fsthalf] = cross_correlation(datameg_xcor_fsthalf, dataphotoL_xcor_fsthalf, dataphotoR_xcor_fsthalf, nTrials/2);
% [xcorL_sndhalf, xcorR_sndhalf, pDelays_sndhalf] = cross_correlation(datameg_xcor_sndhalf, dataphotoL_xcor_sndhalf, dataphotoR_xcor_sndhalf, nTrials/2);
% 
% %plot cross-correlation
% figure;
% subplot(221); plot(pDelays_fsthalf, xcorL_fsthalf); ylim([-.1 .1]); xlabel('delays'); ylabel('correlation');
% subplot(223); plot(pDelays_fsthalf, abs(hilbert(xcorL_fsthalf'))); ylim([0 .1]); xlabel('delays'); ylabel('correlation');
% subplot(222); plot(pDelays_fsthalf, xcorR_fsthalf); ylim([-.1 .1]); xlabel('delays'); ylabel('correlation');
% subplot(224); plot(pDelays_fsthalf, abs(hilbert(xcorR_fsthalf'))); ylim([0 .1]); xlabel('delays'); ylabel('correlation');
% figure;
% subplot(221); plot(pDelays_sndhalf, xcorL_sndhalf); ylim([-.1 .1]); xlabel('delays'); ylabel('correlation');
% subplot(223); plot(pDelays_sndhalf, abs(hilbert(xcorL_sndhalf'))); ylim([0 .1]); xlabel('delays'); ylabel('correlation');
% subplot(222); plot(pDelays_sndhalf, xcorR_sndhalf); ylim([-.1 .1]); xlabel('delays'); ylabel('correlation');
% subplot(224); plot(pDelays_sndhalf, abs(hilbert(xcorR_sndhalf'))); ylim([0 .1]); xlabel('delays'); ylabel('correlation');
% 
% % Topoplot left > right
% %get values
% xcorLhilb_fsthalf = abs(hilbert(xcorL_fsthalf'));
% xcorRhilb_fsthalf = abs(hilbert(xcorR_fsthalf'));
% xcorLhilb_sndhalf = abs(hilbert(xcorL_sndhalf'));
% xcorRhilb_sndhalf = abs(hilbert(xcorR_sndhalf'));
% 
% xcorLhilb_fsthalf_value = xcorLhilb_fsthalf(141,:);
% xcorRhilb_fsthalf_value = xcorRhilb_fsthalf(141,:);
% xcorLhilb_sndhalf_value = xcorLhilb_sndhalf(141,:);
% xcorRhilb_sndhalf_value = xcorRhilb_sndhalf(141,:);
% 
% xcorLmR_fsthalf = (xcorLhilb_fsthalf_value-xcorRhilb_fsthalf_value);%./(xcorLhilb_fsthalf_value+xcorRhilb_fsthalf_value);
% xcorLmR_sndhalf = (xcorLhilb_sndhalf_value-xcorRhilb_sndhalf_value);%./(xcorLhilb_sndhalf_value+xcorRhilb_sndhalf_value);
% 
% %topoplot
% load tfr_all
% tfr_all.freq = 1;
% tfr_all.time = 1;
% tfr_all.dimord = 'chan_freq';
% for i=1:2
%     if i == 1
%         tfr_all.powspctrm = xcorLmR_fsthalf';
%     else
%         tfr_all.powspctrm = xcorLmR_sndhalf';
%     end
%     %combine sensors
%     cfg             = [];
%     cfg.method      = 'sum';
%     tfr_all_c = ft_combineplanar(cfg,tfr_all);  
%     %topoplot
%     cfg = [];
%     cfg.zlim            = [-.15 .15] ;
%     cfg.layout = 'neuromag306cmb.lay';
%     figure;ft_topoplotER(cfg,tfr_all_c)
%     colormap jet
% end
% 
% %topoplot of the difference
% xcorLmR_SmF = xcorLmR_sndhalf - xcorLmR_fsthalf;
% tfr_all.powspctrm = xcorLmR_SmF';
% %combine sensors
% cfg             = [];
% cfg.method      = 'sum';
% tfr_all_c = ft_combineplanar(cfg,tfr_all);  
% %topoplot
% cfg = [];
% cfg.zlim            = [-.02 .02] ;
% cfg.layout = 'neuromag306cmb.lay';
% figure;ft_topoplotER(cfg,tfr_all_c)
% colormap jet


%% Cross vs partial correlation
% diffcorL = xcorL-pcorL;
% diffcorR = xcorR-pcorR;
% %plot
% figure;
% subplot(121); plot(pDelays, diffcorL); ylim([-.2 .2]); xlabel('delays'); ylabel('correlation');
% subplot(122); plot(pDelays, diffcorR); ylim([-.2 .2]); xlabel('delays'); ylabel('correlation');


