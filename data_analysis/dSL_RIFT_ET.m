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


%% Select trials with good fixation

% check trial fixation
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
