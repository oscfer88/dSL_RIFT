% Time-frequency analysis of power

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

%% Reading and segmenting the data creating trials

% Get fequency tagging start trigger
trigger_eventtype = [2,10,18,26,34,42,50,58,66,74,82,90,98,106,114,122,130,138,146,154,162,170,178,186,194,202,210,218,226,234,242,250];

%Define trials and preprocess data
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
    cfg.channel                 = {'MEGGRAD'}; % 'MEGGRAD' 'MISC004' 'MISC005'
    cfg.detrend                 = 'yes';    
    data(i,1) = ft_preprocessing(cfg); 

end

% Merge datafiles
if size(data,1) == 1
    data_planar = ft_appenddata([], data(1,1));
elseif size(data,1) == 2
    data_planar = ft_appenddata([], data(1,1),data(2,1));
elseif size(data,1) == 3
    data_planar = ft_appenddata([], data(1,1),data(2,1),data(3,1));
end
data_planar = rmfield(data_planar,'sampleinfo'); %remove the 'sampleinfo' field of the newly created datafile

%% Artifact rejections in planar gradiometers I
% data_planar.trialinfo(1:end,2)=1:length(data_planar.trialinfo);
% 
% cfg                     = [];
% cfg.channel             = [1:25]; 
% cfg.viewmode            = 'vertical';
% cfg.ylim                = [-4.5e-12 4.5e-12];
% cfg.continuous          = 'no';
% ft_databrowser(cfg, data_planar);
% 
% cfg                     = [];
% cfg.method              = 'summary';
% cfg.layout              = 'neuromag306planar.lay';
% rjv1 = ft_rejectvisual(cfg, data_planar);
% 
% trl_keep = rjv1.trialinfo(1:end,2); %sbj1_s1: 68 trials removed; sbj1_s2: 39 trials removed; sbj2_s1: 43 trials removed; sbj2_s2: 34 trials removed; sbj3_s1: 87 trials removed; sbj3_s2: 116 trials removed; sbj4_s1: 68 trials removed; sbj4_s2: 62 trials removed; sbj5_s1: 163 trials removed; sbj5_s2: 129 trials removed; sbj6_s1: 117 trials removed; sbj6_s2: 39 trials removed; sbj7_s1: 42 trials removed; sbj7_s2: 38 trials removed; sbj8_s1: 168 trials removed; sbj8_s2: 189 trials removed; sbj9_s1: 145 trials removed; sbj9_s2: 150 trials removed; sbj10_s1: 62 trials removed; sbj10_s2: 62 trials removed; sbj11_s1: 26 trials removed; sbj11_s2: 32 trials removed; sbj12_s1: 165 trials removed; sbj12_s2: 187 trials removed; sbj13_s1: 161 trials removed; sbj13_s2: 105 trials removed; sbj14_s1: 70 trials removed; sbj14_s2: 109 trials removed; sbj15_s1: 72 trials removed; sbj15_s2: 59 trials removed; sbj16_s1: 99 trials removed; sbj16_s2: 61 trials removed; sbj17_s1: 54 trials removed; sbj17_s2: 43 trials removed; sbj18_s1: 73 trials removed; sbj18_s2: 71 trials removed; sbj19_s1: 107 trials removed; sbj19_s2: 133 trials removed; sbj20_s1: 36 trials removed; sbj20_s2: 54 trials removed
% rjv1.trialinfo(:,2) = [];
% chan_rej = setdiff(data_planar.label,rjv1.label); %all: none
% chan_keep = data_planar.label;
% chan_keep(find(ismember(data_planar.label(:,1), chan_rej(:,1)))) = [];
% % save('trl_keep','trl_keep')
% % save('chan_keep','chan_keep')

load trl_keep
load good_trials %correct response trials
load chan_keep
load eye_keep2
nolear_trials = 101:720; %remove learning (the first 100 trials)
cfg                     = [];
cfg.trials              = intersect(trl_keep, intersect(eye_keep, intersect(good_trials, nolear_trials)));
cfg.channel             = chan_keep;
planar_rjv = ft_selectdata(cfg,data_planar);

%% Artifact suppression by ICA
% cfg                     = []; 
% cfg.resamplefs          = 300;   
% planar_resamp = ft_resampledata(cfg, planar_rjv); %downsampling
% 
% %step 1 - ICA Decomposition
% cfg                     = [];  
% cfg.method              = 'runica';  
% cfg.runica.maxsteps     = 100;  
% planar_comp  = ft_componentanalysis(cfg, planar_resamp);
% % save('comp','planar_comp')
% 
% %step 2 - identifying components
% load comp3
% cfg = [];
% cfg.channel = [1:10]; 
% cfg.continuous='no';
% cfg.viewmode = 'component'; 
% cfg.layout = 'neuromag306planar.lay';
% ft_databrowser(cfg, planar_comp);

%step 3 - rejecting components
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

planar_clear  = ft_rejectcomponent(cfg, planar_comp, planar_rjv);

% Artifact rejections in planar gradiometers II
% cfg                     = [];
% cfg.channel             = [1:25]; 
% cfg.viewmode            = 'vertical';
% cfg.ylim                = [-4.5e-12 4.5e-12];
% cfg.continuous          = 'no';
% ft_databrowser(cfg, planar_clear);
% 
% cfg                     = [];
% cfg.method              = 'summary';
% cfg.layout              = 'neuromag306planar.lay';
% rjv2 = ft_rejectvisual(cfg, planar_clear);

%% Perform the time-frequency analysis of power; all frequencies
cfg                         = [];
cfg.output                  = 'pow';
cfg.channel                 = 'MEGGRAD';
cfg.taper                   = 'hanning';
cfg.method                  = 'mtmconvol';
cfg.foi                     = 1:1:100;
% numfoi                      = length(cfg.foi);
cfg.t_ftimwin               = ones(length(cfg.foi),1).* 1;
cfg.toi                     = [-1:0.05:2];
cfg.keeptrials              = 'no'; 
tfr_all = ft_freqanalysis(cfg, planar_clear);
% save('tfr_all', 'tfr_all')

%% Plots

%combine the power of two gradiometers
% load tfr_all
cfg             = [];
cfg.method      = 'sum';
tfr_all_c = ft_combineplanar(cfg,tfr_all);

% %baseline correction
% tfr_data_baselined = cell(size(filelist,1),1);
% cfg               = [];
% cfg.baseline      = [-0.5 -0.3];
% cfg.baselinetype  = 'relative';
% for i = 1:size(filelist,1)
%     tfr_data_baselined{i} = ft_freqbaseline(cfg, tfr_data{i});
% end

%plot the TFR results
% cfg                 = []; %plot all planar sensors
% cfg.baseline        = [-0.5 -0.3]; 
% cfg.baselinetype    = 'relative';
% cfg.xlim            = [-0.5 2];
% cfg.zlim            = [0  2] ;	        
% cfg.showlabels      = 'no';	 
% cfg.layout          = 'neuromag306planar.lay';
% ft_multiplotTFR(cfg, tfr_all);
% colormap jet

cfg                 = []; %plot all combined sensors
cfg.baseline        = [-0.5 -0.3]; 
cfg.baselinetype    = 'relative';
cfg.xlim            = [-1 2];
cfg.zlim            = [0  2];	        
cfg.showlabels      = 'no';	 
cfg.layout          = 'neuromag306cmb.lay';
ft_multiplotTFR(cfg, tfr_all_c);
colormap jet

% cfg                 = []; %plot all combined sensors W/O baseline correction
% cfg.xlim            = [-1 2];
% cfg.layout          = 'neuromag306cmb.lay';
% ft_multiplotTFR(cfg, tfr_all_c);
% colormap jet

% cfg                 = []; %plot single sensor
% cfg.baseline        = [-0.5 -0.3]; 
% cfg.baselinetype    = 'relative';
% cfg.xlim            = [-0.5 2];
% cfg.zlim            = [0  2];	        
% cfg.showlabels      = 'no';	
% cfg.channel         = 'MEG2112+2113'; 
% ft_singleplotTFR(cfg, tfr_all_c);
% colormap jet

cfg                 = []; %plot topography of alpha (8-12 HZ)
cfg.baseline        = [-0.5 -0.3]; 
cfg.baselinetype    = 'relative';
cfg.xlim            = [0 1.5];
cfg.ylim            = [8 12];
cfg.zlim            = [0  2];	        
cfg.showlabels       = 'no';	
cfg.layout          = 'neuromag306cmb.lay';
ft_topoplotTFR(cfg, tfr_all_c);
colormap jet

% cfg                 = []; %plot topography of alpha (8-12 HZ) as a function of time
% cfg.baseline        = [-0.5 -0.3]; 
% cfg.baselinetype    = 'relative';
% cfg.xlim            = [-1:0.5:2];
% cfg.ylim            = [8 12];
% cfg.zlim            = [0  2];	        
% cfg.showlabels       = 'no';	
% cfg.layout          = 'neuromag306cmb.lay';
% ft_topoplotTFR(cfg, tfr_all_c);
% colormap jet

% cfg                 = []; %plot topography of alpha (8-12 HZ) as a function of time W/O baseline correction
% cfg.xlim            = [-1:0.5:2];
% cfg.ylim            = [8 12];
% cfg.layout          = 'neuromag306cmb.lay';
% ft_topoplotTFR(cfg, tfr_all_c);
% colormap jet

cfg                 = []; %plot topography of bbt (50-80 HZ)
cfg.baseline        = [-0.5 -0.3]; 
cfg.baselinetype    = 'relative';
cfg.xlim            = [0 1.5];
cfg.ylim            = [50 80];
cfg.zlim            = [0.9  1.1];	        
cfg.showlabels       = 'no';	
cfg.layout          = 'neuromag306cmb.lay';
ft_topoplotTFR(cfg, tfr_all_c);
colormap jet

% cfg                 = []; %plot topography of bbt (50-80 HZ) as a function of time
% cfg.baseline        = [-0.5 -0.3]; 
% cfg.baselinetype    = 'relative';
% cfg.xlim            = [-0.4:0.4:2];
% cfg.ylim            = [50 80];
% cfg.zlim            = [0  2];	        
% cfg.showlabels       = 'no';	
% cfg.layout          = 'neuromag306cmb.lay';
% ft_topoplotTFR(cfg, tfr_all_c);
% colormap jet


%%%%%%%%%%%%%%%%%%%%%%
% % Divide trials in two blocks
% idx = round(quantile(1:size(planar_clear.trial,2),[0 0.50 1.00]));
% planar_block = cell(2,1);
% 
% cfg = [];
% cfg.trials = (idx(1):idx(2))';
% planar_block{1} = ft_selectdata(cfg, planar_clear);
% 
% cfg = [];
% cfg.trials = (1+idx(1):idx(2))';
% planar_block{2} = ft_selectdata(cfg, planar_clear);
% 
% cfg = [];
% cfg.trials = (1+idx(2):idx(3))';
% planar_block{3} = ft_selectdata(cfg, planar_clear);
% 
% cfg = [];
% cfg.trials = (1+idx(3):idx(4))';
% planar_block{4} = ft_selectdata(cfg, planar_clear);
% 
% % Perform the time-frequency analysis of power; all frequencies
% tfr_block = cell(4,1);
% for i= 1:4
% 	cfg                         = [];
% 	cfg.output                  = 'pow';
%     cfg.channel                 = 'MEGGRAD';
%     cfg.taper                   = 'hanning';
%     cfg.method                  = 'mtmconvol';
%     cfg.foi                     = 1:1:100;
%     numfoi                      = length(cfg.foi);
%     cfg.t_ftimwin               = ones(length(cfg.foi),1).* 1;
%     cfg.toi                     = [-1.5 : 0.05: 2];
%     cfg.keeptrials              = 'no'; 
%     tfr_block{i} = ft_freqanalysis(cfg, planar_block{i}); 
% end
% % save('tfr_block','tfr_block');
% 
% %combine the power of two gradiometers
% tfr_block_c = cell(4,1);
% for i=1:4
%     cfg             = [];
%     cfg.method      = 'sum';
%     tfr_block_c{i} = ft_combineplanar(cfg,tfr_block{i});
% end
% 
% %plot the TFR results
% for i=1:4
%     figure;
%     cfg                 = []; %plot all combined sensors
%     cfg.baseline        = [-0.5 -0.3]; 
%     cfg.baselinetype    = 'relative';
%     cfg.xlim            = [-0.5 2];
%     cfg.zlim            = [0.5  1.5] ;	        
%     cfg.layout          = 'neuromag306cmb.lay';
%     ft_multiplotTFR(cfg, tfr_block_c{i});
% end
% 
% for i=1:4
%     figure;
%     cfg                 = []; %plot topography of alpha (8-12 HZ) as a function of time
%     cfg.baseline        = [-0.5 -0.3]; 
%     cfg.baselinetype    = 'relative';
%     cfg.xlim            = [-0.5:0.5:2];
%     cfg.ylim            = [8 12];
%     cfg.zlim            = [0.5  1.5] ;	        
%     cfg.layout          = 'neuromag306cmb.lay';
%     ft_topoplotTFR(cfg, tfr_block_c{i});
% end
% 
% for i=1:4
%     figure;
%     cfg                 = []; %plot topography of bbt (55-75 HZ) as a function of time
%     cfg.baseline        = [-0.5 -0.3]; 
%     cfg.baselinetype    = 'relative';
%     cfg.xlim            = [-0.5:0.5:2];
%     cfg.ylim            = [55 75];
%     cfg.zlim            = [0.5  1.5];	        
%     cfg.layout          = 'neuromag306cmb.lay';
%     ft_topoplotTFR(cfg, tfr_block_c{i});
% end
% 
% 
% % Compute difference between first block and other blocks
% tfr_block_c_diff = cell(3,1);
% for i=1:3
%     tfr_block_c_diff{i} = tfr_block_c{1};
%     tfr_block_c_diff{i}.powspctrm = (tfr_block_c{i+1}.powspctrm-tfr_block_c{1}.powspctrm)./(tfr_block_c{i+1}.powspctrm+tfr_block_c{1}.powspctrm);
% end
% 
% %plot the TFR results
% for i=1:3
%     figure;
%     cfg                 = []; %plot all combined sensors
%     cfg.xlim            = [-0.5 2];
%     cfg.zlim            = [-0.15 0.15] ;	        
%     cfg.layout          = 'neuromag306cmb.lay';
%     ft_multiplotTFR(cfg, tfr_block_c_diff{i});
% end
% 
% for i=1:3
%     cfg                 = []; %plot topography of alpha (8-12 HZ) as a function of time
%     cfg.xlim            = [0 1.5];
%     cfg.ylim            = [8 12];
%     cfg.zlim            = [-0.15 0.15] ;	        
%     cfg.layout          = 'neuromag306cmb.lay';
%     subplot(1,3,i)
%     ft_topoplotTFR(cfg, tfr_block_c_diff{i});
% end
% 
% for i=1:3
%     cfg                 = []; %plot topography of bbt (55-75 HZ) as a function of time
%     cfg.xlim            = [0 1.5];
%     cfg.ylim            = [55 75];
%     cfg.zlim            = [-0.15 0.15];	        
%     cfg.layout          = 'neuromag306cmb.lay';
%     subplot(1,3,i)
%     ft_topoplotTFR(cfg, tfr_block_c{i});
% end
% %%%%%%%%%%%%%%%%%%%%%%





% % Left vs right sensors
% cfg                     = [];
% cfg.channel             = {
%     'MEG0722' 'MEG0723' 'MEG0732' 'MEG0733' ...
%     'MEG0912' 'MEG0913' 'MEG0922' 'MEG0923' 'MEG0932' 'MEG0933' 'MEG0942' 'MEG0943' ...
%     'MEG1022' 'MEG1023' 'MEG1032' 'MEG1033' 'MEG1042' 'MEG1043' ...
%     'MEG1112' 'MEG1113' 'MEG1122' 'MEG1123' 'MEG1132' 'MEG1133' 'MEG1142' 'MEG1143' ...
%     'MEG1212' 'MEG1213' 'MEG1222' 'MEG1223' 'MEG1232' 'MEG1233' 'MEG1242' 'MEG1243' ...
%     'MEG1312' 'MEG1313' 'MEG1322' 'MEG1323' 'MEG1332' 'MEG1333' 'MEG1342' 'MEG1343' ...
%     'MEG1412' 'MEG1413' 'MEG1422' 'MEG1423' 'MEG1432' 'MEG1433' 'MEG1442' 'MEG1443' ...
%     'MEG2022' 'MEG2023' 'MEG2032' 'MEG2033' ...
%     'MEG2132' 'MEG2133' ...
%     'MEG2212' 'MEG2213' 'MEG2222' 'MEG2223' 'MEG2232' 'MEG2233' 'MEG2242' 'MEG2243' ...
%     'MEG2312' 'MEG2313' 'MEG2322' 'MEG2323' 'MEG2332' 'MEG2333' 'MEG2342' 'MEG2343' ...
%     'MEG2412' 'MEG2413' 'MEG2422' 'MEG2423' 'MEG2432' 'MEG2433' 'MEG2442' 'MEG2443' ...
%     'MEG2512' 'MEG2513' 'MEG2522' 'MEG2523' 'MEG2532' 'MEG2533' 'MEG2542' 'MEG2543' ...
%     'MEG2612' 'MEG2613' 'MEG2622' 'MEG2623' 'MEG2632' 'MEG2633' 'MEG2642' 'MEG2643'};
% tfr_right = ft_selectdata(cfg, tfr_all);
% 
% cfg                     = [];
% cfg.channel             = {
%     'MEG0112' 'MEG0113' 'MEG0122' 'MEG0123' 'MEG0132' 'MEG0133' 'MEG0142' 'MEG0143' ...
%     'MEG0212' 'MEG0213' 'MEG0222' 'MEG0223' 'MEG0232' 'MEG0233' 'MEG0242' 'MEG0243' ...
%     'MEG0312' 'MEG0313' 'MEG0322' 'MEG0323' 'MEG0332' 'MEG0333' 'MEG0342' 'MEG0343' ...
%     'MEG0412' 'MEG0413' 'MEG0422' 'MEG0423' 'MEG0432' 'MEG0433' 'MEG0442' 'MEG0443' ...
%     'MEG0512' 'MEG0513' 'MEG0522' 'MEG0523' 'MEG0532' 'MEG0533' 'MEG0542' 'MEG0543' ...
%     'MEG0612' 'MEG0613' 'MEG0632' 'MEG0633' 'MEG0642' 'MEG0643' ...
%     'MEG0712' 'MEG0713' 'MEG0742' 'MEG0743' ...
%     'MEG1512' 'MEG1513' 'MEG1522' 'MEG1523' 'MEG1532' 'MEG1533' 'MEG1542' 'MEG1543' ...
%     'MEG1612' 'MEG1613' 'MEG1622' 'MEG1623' 'MEG1632' 'MEG1633' 'MEG1642' 'MEG1643' ...
%     'MEG1712' 'MEG1713' 'MEG1722' 'MEG1723' 'MEG1732' 'MEG1733' 'MEG1742' 'MEG1743' ...
%     'MEG1812' 'MEG1813' 'MEG1822' 'MEG1823' 'MEG1832' 'MEG1833' 'MEG1842' 'MEG1843' ...
%     'MEG1912' 'MEG1913' 'MEG1922' 'MEG1923' 'MEG1932' 'MEG1933' 'MEG1942' 'MEG1943' ...
%     'MEG2012' 'MEG2013' 'MEG2042' 'MEG2043' ...
%     'MEG2142' 'MEG2143'};
% tfr_left = ft_selectdata(cfg, tfr_all);
