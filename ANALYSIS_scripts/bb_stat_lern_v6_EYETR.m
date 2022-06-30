%% ASCII

clearvars

% Load data
folder = 'Q:\EYELINK_data';
files = dir([folder, '\*.asc']);

% Get trial info
trials = cell(size(files,1),1);
for i=1:size(files,1)
    
    fprintf(1,'#%d\n',i); %print file nr

    % load data
    data = read_eyelink_asc([folder,'\',files(i).name]);
    
    % define trials
    trials{i}.mess = data.msg(contains(data.msg, 'Distractor'));
    %trials{i}.idx = find(contains(data.msg, 'Distractor'));
    trials{i}.time = str2double(extractBetween(trials{i}.mess, 'MSG	', ' '));
    for j=1:size(trials{i}.time,1)
        trials{i}.dottoi(j,:) = [trials{i}.time(j), trials{i}.time(j)+499]; %start and end of time of fixation dot screen
        IDXdotpos = find(data.dat(1,:) == trials{i}.dottoi(j,1)):find(data.dat(1,:) == trials{i}.dottoi(j,2));
        trials{i}.dotxpos(j,:) = data.dat(2,IDXdotpos);
        trials{i}.dotypos(j,:) = data.dat(3,IDXdotpos);
        trials{i}.tagtoi(j,:) = [trials{i}.time(j)+500, trials{i}.time(j)+1999]; %start and end of pre-search array
    end
end
% save('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data\Eyelink\trials', 'trials')

% Get event info
events = cell(size(files,1),1);
for i=1:size(files,1)
    % define events
    %fixations
    events{i}.fixmess = data.efix;
    split_mess = [];
    for j=1:size(events{i}.fixmess,1)
        split_mess = strsplit(events{i}.fixmess{j});
        events{i}.fixtoi(j,:) = [str2double(split_mess(3)), str2double(split_mess(4))];
        events{i}.fixdur(j,:) = str2double(split_mess(5));
        events{i}.fixapos(j,:) = [str2double(split_mess(6)), str2double(split_mess(7))]; %averaged xy position
    end
    %blinks
    events{i}.blinkmess = data.eblink;
    split_mess = [];
    for j=1:size(events{i}.blinkmess,1)
        split_mess = strsplit(events{i}.blinkmess{j});
        events{i}.blinktoi(j,:) = [str2double(split_mess(3)), str2double(split_mess(4))];
        events{i}.blinkdur(j,:) = str2double(split_mess(5));
    end
    %saccades
    events{i}.sacmess = data.esacc;
    split_mess = [];
    events{i}.sactoi = [];
    for j=1:size(events{i}.sacmess,1)
        split_mess = strsplit(events{i}.sacmess{j});
        events{i}.sactoi(j,:) = [str2double(split_mess(3)), str2double(split_mess(4))]; %stime etime deg
        events{i}.sacdur(j,:) = str2double(split_mess(5));
        events{i}.sacspos(j,:) = [str2double(split_mess(6)), str2double(split_mess(7))];
        events{i}.sacepos(j,:) = [str2double(split_mess(8)), str2double(split_mess(9))];
        events{i}.sacampl(j,:) = str2double(split_mess(10)); %in dergees
    end
    events{i}.sacmess = data.esacc;
    
end
% save('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data\Eyelink\events', 'events')

% Divide events into trials
events_by_trials = cell(size(trials,1),1);
for i=1:size(trials,1) %loop over subject files
    
    % drift correction
    DOT = trials{i}.dottoi;
    DOTX = trials{i}.dotxpos;
    DOTY = trials{i}.dotypos;
    BLINK = events{i}.blinktoi;
    for j=1:size(DOT,1)
        DOTBLINK = zeros(1,size(DOT(1,1):DOT(1,2),2));
        for k=1:size(BLINK,1)
            DOTBLINK0 = ismember(DOT(j,1):DOT(j,2), BLINK(k,1):BLINK(k,2));
            DOTBLINK = DOTBLINK | DOTBLINK0;
        end
        DOTX(j,DOTBLINK) = nan;
        DOTY(j,DOTBLINK) = nan;
    end
    %get median fixation position and drift
    DRIFT = [nanmedian(DOTX,2), nanmedian(DOTY,2)];
    DRIFTCOR = [scrW/2, scrH/2] - DRIFT;
    if isnan(DRIFTCOR(1,:))
        DRIFTCOR(1,:) = [0 0];
    end
    while any(isnan(DRIFTCOR(:,1)))
        DRIFTCORnan = find(isnan(DRIFTCOR(:,1)));
        DRIFTCOR(DRIFTCORnan,:) = DRIFTCOR(DRIFTCORnan-1,:);
    end
    
    % save variables
    events_by_trials{i}.drift = DRIFTCOR';
    events_by_trials{i}.dotxpos = (trials{i}.dotxpos + DRIFTCOR(:,1))';
    events_by_trials{i}.dotypos = (trials{i}.dotypos + DRIFTCOR(:,2))';
    
    for j=1:size(trials{i}.mess,1) %loop over trials
        
        TRIAL = trials{i}.tagtoi(j,1):trials{i}.tagtoi(j,2);
        
        FIXIDX = ismember(events{i}.fixtoi, TRIAL);
        FIXIDX = any(FIXIDX,2);
        events_by_trials{i}.fixtoi{j} = events{i}.fixtoi(FIXIDX,:);
        
        events_by_trials{i}.fixdur{j} = events{i}.fixdur(FIXIDX);
        events_by_trials{i}.fixapos{j} = events{i}.fixapos(FIXIDX,:);
        
        BLINKIDX = ismember(events{i}.blinktoi, TRIAL);
        BLINKIDX = any(BLINKIDX,2);
        events_by_trials{i}.blinktoi{j} = events{i}.blinktoi(BLINKIDX,:);
        events_by_trials{i}.blinkdur{j} = events{i}.blinkdur(BLINKIDX);
        
        if ~isempty(events{i}.sacmess)
            SACIDX = ismember(events{i}.sactoi, TRIAL);
            SACIDX = any(SACIDX,2);
            events_by_trials{i}.sactoi{j} = events{i}.sactoi(SACIDX,:);
            events_by_trials{i}.sacdur{j} = events{i}.sacdur(SACIDX);
            events_by_trials{i}.sacspos{j} = events{i}.sacspos(SACIDX,:);
            events_by_trials{i}.sacepos{j} = events{i}.sacepos(SACIDX,:);
            events_by_trials{i}.sacampl{j} = events{i}.sacampl(SACIDX);
        end
        
    end
    
end
% save('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data\Eyelink\events_by_trials', 'events_by_trials')


%% Saccade analysis

% plot saccades
load('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data\Eyelink\trials')
load('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data\Eyelink\events')

scrW = 1920;
scrWcm = 70.6;
scrH = 1080;
scrDEGcm = 2.5;
scrDEG = scrW/scrWcm*scrDEGcm;

selsacc = cell(size(trials,1),1);
for i=1:size(trials,1)
    
    fprintf(1,'#%d\n',i); %print file nr
    
    if ~isempty(events{i}.sacmess)
        
        % get index of saccades during trial time window
        TAG = trials{i}.tagtoi;
        KEEPTAG = [];
        for j=1:size(events{i}.sactoi,1)
            for k=1:size(TAG,1)
                if events{i}.sactoi(j,1) > TAG(k,1) && events{i}.sactoi(j,2) < TAG(k,2)
                    KEEPTAG = [KEEPTAG, j];
                end
            end
        end
        
        % get index of saccades with blinks in between and before or after
        % 100 ms from a blink
        BLINK = events{i}.blinktoi;
        REMBLINK = [];
        for j=1:size(events{i}.sactoi,1)
            for k=1:size(BLINK,1)
                if ( events{i}.sactoi(j,1) < BLINK(k,1) && events{i}.sactoi(j,2) > BLINK(k,2) ) || ...
                        ( (BLINK(k,2) - events{i}.sactoi(j,1)) > 0 && (BLINK(k,2) - events{i}.sactoi(j,1)) < 100 ) || ...
                        ( (events{i}.sactoi(j,2) - BLINK(k,1)) > 0 && (events{i}.sactoi(j,2) - BLINK(k,1)) < 100 )
                    REMBLINK = [REMBLINK, j];
                end
            end
        end
        
        % get index of very quick and short (micro)saccades
        MICROS = find(events{i}.sacdur < 16 | events{i}.sacampl < 2);
        
        % saccade info and reject trials
        IDX = setdiff(setdiff(KEEPTAG, REMBLINK), MICROS);
        
        SS = events{i}.sacspos(IDX,:); %saccade start position
        SS2 = [];
        SE = events{i}.sacepos(IDX,:); %saccade end position
        SE2 = [];
        STIME = events{i}.sactoi(IDX,1);
        ETIME = events{i}.sactoi(IDX,2);
        DUR = events{i}.sacdur(IDX);
        SIZE = events{i}.sacampl(IDX);
        
        if ~isempty(SS)
            
            % apply drift correction
            DOT = trials{i}.dottoi;
            DOTX = trials{i}.dotxpos;
            DOTY = trials{i}.dotypos;
            %remove blinks from fixation
            DOTXFIX = DOTX;
            DOTYFIX = DOTY;
            for j=1:size(DOT,1)
                DOTBLINK = zeros(1,size(DOT(1,1):DOT(1,2),2));
                for k=1:size(BLINK,1)
                    DOTBLINK0 = ismember(DOT(j,1):DOT(j,2), BLINK(k,1):BLINK(k,2));
                    DOTBLINK = DOTBLINK | DOTBLINK0;
                end
                DOTXFIX(j,DOTBLINK) = nan;
                DOTYFIX(j,DOTBLINK) = nan;
            end
            %get median fixation position and drift
            DRIFT = [nanmedian(DOTXFIX,2), nanmedian(DOTYFIX,2)];
            DRIFTCOR = [scrW/2, scrH/2] - DRIFT;
            if isnan(DRIFTCOR(1,:))
                DRIFTCOR(1,:) = [0 0];
            end
            while any(isnan(DRIFTCOR(:,1)))
                DRIFTCORnan = find(isnan(DRIFTCOR(:,1)));
                DRIFTCOR(DRIFTCORnan,:) = DRIFTCOR(DRIFTCORnan-1,:);
            end
            %actual drift correction
            for j=1:size(TAG,1)
                for k=1:size(SS,1)
                    if STIME(k) > TAG(j,1) && ETIME(k) < TAG(j,2) %check time
                        SS2(k,1) = SS(k,1) + DRIFTCOR(j,1);
                        SS2(k,2) = SS(k,2) + DRIFTCOR(j,2);
                        SE2(k,1) = SE(k,1) + DRIFTCOR(j,1);
                        SE2(k,2) = SE(k,2) + DRIFTCOR(j,2);
                    end
                end
            end
            
            % remove saccade falling out of the screen
            OUTSX = SS2(:,1) < 0 | SS2(:,1) > scrW;
            OUTSY = SS2(:,2) < 0 | SS2(:,2) > scrH;
            OUTEX = SE2(:,1) < 0 | SE2(:,1) > scrW;
            OUTEY = SE2(:,2) < 0 | SE2(:,2) > scrH;
            OUT = ~(OUTSX | OUTSY | OUTEX | OUTEY);
            
            SS2 = SS2(OUT,:);
            SE2 = SE2(OUT,:);
            STIME = STIME(OUT);
            ETIME = ETIME(OUT);
            DUR = DUR(OUT);
            SIZE = SIZE(OUT);
            
            if ~isempty(SS2)

                % write variable
                selsacc{i}.pos = [SS2(:,1), SE2(:,1), SS2(:,2), SE2(:,2)];
                selsacc{i}.time = [STIME, ETIME];
                selsacc{i}.dur = DUR;
                selsacc{i}.size = SIZE;

%                 % plot
%                 for j=1:size(SS2,1)
%                     plot([SS2(j,1),SE2(j,1)],[SS2(j,2),SE2(j,2)]); hold on
%                 end
%                 xlim([0 1920]); ylim([0 1080]); title(['#', num2str(i)]); set(gca,'ydir','reverse'); %used when plotting
%                 line([scrW/2-scrDEG*3, scrW/2+scrDEG*3, scrW/2+scrDEG*3, scrW/2-scrDEG*3, scrW/2-scrDEG*3], ...
%                     [scrH/2-scrDEG*3, scrH/2-scrDEG*3, scrH/2+scrDEG*3, scrH/2+scrDEG*3, scrH/2-scrDEG*3]);
%                 hold off
%                 pause %used when plotting
                
            end
            
        end
                
    end
    
end
% save('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data\Eyelink\selsacc', 'selsacc')

% divide by subject
load('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data\Eyelink\selsacc')
sbj_list = repmat(1:40,2,1);
sbj_list = sbj_list(:);
sbj_list(4) = [];
selsacc_sbj = cell(40,1);
sbj = [];
for i=1:size(selsacc,1)
    if sbj_list(i) == sbj
        temp2 = [temp; selsacc{i}];
        temp = temp2;
        selsacc_sbj{sbj} = [];
    else
        temp = selsacc{i};
    end
    sbj = sbj_list(i);
    selsacc_sbj{sbj} = temp;
end

% get saccades polar coordinates
for i=1:size(selsacc_sbj,1)
    
    fprintf(1,'#%d\n',i); %print file nr
    
    if ~isempty(selsacc_sbj{i})
        
        DS = selsacc_sbj{i}.pos;
        
        XC = DS(:,1);
        YC = DS(:,3);
        
        DS(:,1:2) = DS(:,1:2) - XC;
        DS(:,3:4) = DS(:,3:4) - YC;
        
        %transform to polar
        XS = DS(:,1:2);
        YS = DS(:,3:4);
        [THETA,RHO] = cart2pol(XS, YS); %angular and radial coords
        selsacc_sbj{i,2} = DS;
        selsacc_sbj{i,3} = [THETA, RHO];
        
%         %plot saccade hist
%         polarhistogram(THETA(:,2),180);
%         title(['#', num2str(i)]); %rlim([0 30]);
%         pause
        
    end
    
end


% Divide session by distractor probability lateralization
selsacc_sbjHiLo = selsacc_sbj([1 4 5 8  9 12 13 16 17 20 21 24 25 28 29 32 33 36 37 40],:); %HIGH-LOW condition
selsacc_sbjLoHi = selsacc_sbj([2 3 6 7 10 11 14 15 18 19 22 23 26 27 30 31 34 35 38 39],:); %LOW-HIGH condition

% Divide saccades into leftwards and rightwards
for i=1:size(selsacc_sbjHiLo,1)
    if ~isempty(selsacc_sbjHiLo{i,3})
        RADANGLE = selsacc_sbjHiLo{i,3}(:,2);
        HiLo_LEFT(i) = numel(RADANGLE(RADANGLE > -pi/2 & RADANGLE < pi/2));
        HiLo_RIGHT(i) = numel(RADANGLE(RADANGLE < -pi/2 | RADANGLE > pi/2));
    end
    if ~isempty(selsacc_sbjLoHi{i,3})
        RADANGLE = selsacc_sbjLoHi{i,3}(:,2);
        LoHi_LEFT(i) = numel(RADANGLE(RADANGLE > -pi/2 & RADANGLE < pi/2));
        LoHi_RIGHT(i) = numel(RADANGLE(RADANGLE < -pi/2 | RADANGLE > pi/2));
    end
end

% [~,p] = ttest(HiLo_LEFT, HiLo_RIGHT) %t-test
% [~,p] = ttest(LoHi_LEFT, LoHi_RIGHT) %t-test



% Kolmogorov-Smirnov test of saccade direction (w/ histogram)
HiLo_hist = [];
LoHi_hist = [];
for i=1:size(selsacc_sbjHiLo,1)
    if ~isempty(selsacc_sbjHiLo{i,3})
        HiLo_hist = [HiLo_hist; selsacc_sbjHiLo{i,3}(:,2)];
    end
    if ~isempty(selsacc_sbjLoHi{i,3})
        LoHi_hist = [LoHi_hist; selsacc_sbjLoHi{i,3}(:,2)];
    end
end

[h,p,ks2stat] = kstest2(HiLo_hist, LoHi_hist) %Kolmogorov-Smirnov test

%plot
figure; polarhistogram(HiLo_hist,180);
figure; polarhistogram(LoHi_hist,180);

% Get saccade laterality by probability condition
HiLo = (HiLo_LEFT - HiLo_RIGHT) ./ (HiLo_LEFT + HiLo_RIGHT);
LoHi = (LoHi_LEFT - LoHi_RIGHT) ./ (LoHi_LEFT + LoHi_RIGHT);

[~,p,~,stats] = ttest(HiLo, LoHi) %t-test

%plot
figure;hold on; plot([HiLo; LoHi],'-.k');
plot([nanmean(HiLo), nanmean(LoHi)]','r','LineWidth',2); hold off;
% boxplot([squeeze(xcorH_sbj), squeeze(xcorL_sbj)]);
xlim([0.5,2.5]); xticks([1,2]); xticklabels({'High-low','Low-high'}); ylim([-1.5 1.5]); ylabel('Saccade direction');

%% Create eye_keep files for trial exclusion V1

% check trial fixation
load('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data\Eyelink\events_by_trials')

scrW = 1920;
scrWcm = 70.6;
scrH = 1080;
scrDEGcm = 2.5;
scrDEG = scrW/scrWcm*scrDEGcm;
eyekeep = nan(size(events_by_trials,1), size(events_by_trials{i}.fixtoi,2));
for i=1:size(events_by_trials,1) %loop over subject files
    
    fprintf(1,'#%d\n',i); %print subject file nr
    
    for j=1:size(events_by_trials{i}.fixtoi,2) %loop over trials
        
        DRIFT = events_by_trials{i}.drift(:,j)';
        FIX = events_by_trials{i}.fixapos{j};
        FIX = FIX + DRIFT;
        
        if isempty(FIX) ||...
                any(FIX(:,1) > (scrW/2 - scrDEG*2)) && any(FIX(:,1) < (scrW/2 + scrDEG*2)) &&...
                any(FIX(:,2) > (scrH/2 - scrDEG*2)) && any(FIX(:,2) < (scrH/2 + scrDEG*2))
            eyekeep(i,j) = 1;
        else
            eyekeep(i,j) = 0;
        end
        
    end
    
end

% get keet trial list
sbj_list = repmat(1:40,2,1);
sbj_list = sbj_list(:);
sbj_list(4) = [];
fix_trials = cell(40,1);
sbj = [];
for i=1:size(eyekeep,1)
    if sbj_list(i) == sbj
        eye2 = [eye, find(eyekeep(i,:)) + 360];
        eye = eye2;
        fix_trials{sbj} = [];
    else
        eye = find(eyekeep(i,:));
    end
    sbj = sbj_list(i);
    fix_trials{sbj} = eye;
end

% save subject X session datafiles
folderMEG = 'Q:/MEG_data/';
filesMEG = dir(fullfile(folderMEG, '**', 'trl_keep.mat'));
filesMEG_ordered = filesMEG([1 7 2 5 3 6 4 8 9 11 10 13 12 15 14 16 17 19 18 20 21 23 22 24 25 28 26 29 27 31 30 34 32 35 33 36 37 39 38 40]);
for i=1:size(filesMEG_ordered,1)
    eye_keep = fix_trials{i};
    save([filesMEG_ordered(i).folder, '\eye_keep2'], 'eye_keep')
end

% %% Create eye_keep files for trial exclusion V0
% 
% % check trial fixation
% load('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data\Eyelink\trials')
% load('C:\Users\ferranto\DataShare\UoB\Projects\01 SL & frequency tagging\Data\Eyelink\events')
% for i=1:size(trials,1)
%     
%     fprintf(1,'#%d\n',i); %print file nr
%     
% %     %V1 - using fixation period
% %     fix = [];
% %     for j=1:size(events{i}.fixtoi,1)
% %         fix = [fix, events{i}.fixtoi(j,1):events{i}.fixtoi(j,2)];
% %     end
% %     for j=1:size(trials{i}.tagtoi,1)    
% %         if mean(ismember([trials{i}.tagtoi(j,1):trials{i}.tagtoi(j,2)],fix))
% %             trials{i}.val(j) = 1;
% %         else
% %             trials{i}.val(j) = 0;
% %         end
% %     end
%     
%     %V2 - using saccade period
%     sac = [];
%     for j=1:size(events{i}.sactoi,1)
%         if events{i}.sactoi(j,3) > 10 && events{i}.sactoi(j,3) < 100 && events{i}.sactoi(j,4) > 2 %<- 3 deg in eye_keep2, used in comp3
%             sac = [sac, events{i}.sactoi(j,1):events{i}.sactoi(j,2)];
%         end
%     end
%     for j=1:size(trials{i}.tagtoi,1)    
%         if any(ismember([trials{i}.tagtoi(j,1):trials{i}.tagtoi(j,2)],sac))
%             trials{i}.val(j) = 0;
%         else
%             trials{i}.val(j) = 1;
%         end
%     end
%     
% end
% 
% % get kept trials
% sbj_list = repmat(1:40,2,1);
% sbj_list = sbj_list(:);
% sbj_list(4) = [];
% fix_trials = cell(40,1);
% sbj = [];
% for i=1:size(trials,1)
%     if sbj_list(i) == sbj
%         eye2 = [eye, find(trials{i}.val) + 360];
%         eye = eye2;
%         fix_trials{sbj} = [];
%     else
%         eye = find(trials{i}.val);
%     end
%     sbj = sbj_list(i);
%     fix_trials{sbj} = eye;
% end
% 
% % save subject X session datafiles
% folderMEG = 'Q:/MEG_data/';
% filesMEG = dir(fullfile(folderMEG, '**', 'trl_keep.mat'));
% filesMEG_ordered = filesMEG([1 7 2 5 3 6 4 8 9 11 10 13 12 15 14 16 17 19 18 20 21 23 22 24 25 28 26 29 27 31 30 34 32 35 33 36 37 39 38 40]);
% for i=1:size(filesMEG_ordered,1)
%     eye_keep = fix_trials{i};
%     save([filesMEG_ordered(i).folder, '\eye_keep'], 'eye_keep')
% end
