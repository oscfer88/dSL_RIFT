%% BBT_SL
%
% Full title: BroadBand frequency Tagging & Statistical Learning of distractor location
% Author: Oscar Ferrante
% 
% 
% PARADIGM:
% This experiment consists in an additional singleton visual search task.
% Four (4) gabor patches are presented across four (4) locations. After a
% short presentation of identical vertically oriented gabors, one patch
% is tilted of 15 degrees leftwards or rightwards. Participants have to
% indicated the direction of this tilt. On 66% of the trials, a distractor
% stimulus is presented in place of one of the remaining vertical patches.
% 'B/W' VERSION: 
% The distractor is a horizontally oriented gabor.
% 'COLOR' VERSION:io
% The distractor is a gabor presented in a different
% color than the taret (target green, distractor red, and vice versa).
% 
% TASK:
% Indicate the orientation of the target while ignoring the distractor.
% 
% SL MANIPULATION:
% When presented, the distractor is displayed in one hemifield on 75% of 
% the trials and in the other on 25% (counterbalanced across participants).
% 
% LOCATIONS:
%               1       3
%                   +
%               2       4
%
% HEMIFIELD COUNTERBALANCE:
% The association between the spatial probability manipulation and the
% hemifield is counterbalanced across participant and session.
% - ODD participants (e.g., 1, 3, etc): the distractor is presented with
% high probability on the LEFT hemifield (loc. 1 and 2) in session 1 and
% viceversa in session 2.
% - EVEN participants (e.g., 2, 4, etc) have the high probability location
% associated with the RIGHT hemifield (loc. 3 and 4) in session 1 and
% viceversa in session 2.
% 
% BROADBAND FREQUENCY TAGGING
% The two stimuli on the bottom (loc 2 and 4) are flickered using two
% different sinusoidal signals composed by frequencies between 45/50 and 85/80 Hz.
% 
% TRIGGER CODE & MISC:
% STI001 - trial onset 
% STI002 - frequency tagging onset and trial info (STI004, STI005, STI006, STI007, STI008)
% STI003 - search array onset
% MISC004 = SENSOR_1 = RIGHT
% MISC005 = SENSOR_2 = LEFT
% 
% 
% CHANGELOG:
% v10:
%     - new broadband frequency tagging signal generation (thanks Alex)
%     - inserted exp_part variable in the dialog box
%     - interted tag_left and tag_right variables in the text log file
%     - moved dur_tag_propixx inside bbTaggingSignal subfunction
%     - new trigger coding:
%           - used a new coding system for distractor location
%           - moved trial info in the frequency tagging onset trigger
%           - removed tag info trigger
%     - new fixation point
%     - added practice epoch and moved practice variable in the dialog box
%     - changed fixation screen duration from 300 ms to 500 ms
%     - changed ITI jitter deviation from 400 ms to 500 ms
%     - blended the border of b\w gabor patches
%     - used the same background color for both versions (color was 0.5)
%     - changed white color for 1 to 0.8
% v11:
%     - corrected the range of the frequency tagging from 45-65 to 55-75 Hz
%     - simplified the gabor function
%     - new color stimuli (color/white gabors)
%     - inserted exp_part as part of the output file name
%     - changed the color of the photodiode pads from background to black
%     - reduced the size of the photodiode pads from 2deg to 50px
% v13:
%     - new color stimuli (light-color/black gabors)
%     - corrected eyelink output filename (now including exp_part)
%     - inserted subject_nr variable
%     - subject_id is now the MEG identification number
%     - inserted exp_part variable in the log file
% v14:
%     - fixed existing output file check [sessionInfo]
%     - corrected all output file names
% v15:
%     - tilted the distractor also in the color version
%     - corrected subject_nr usage
%     - new broadband frequency tagging function [bbTaggingSignal]
%         - the function is now used off-line
%         - a different couple of signals is now produced for each trial
%         - the central frequency f0 is now varied in the range 65-69Hz
%         - inserted a two-step verification during the signal construction
%             - checked if auto-correlations are similar to a noise signal
%             - checed if cross-correlation between signals is lower than 0.1
%     - introduced exp_session variable [sessionInfo]
%         - swaped SL manipulation on second session [trialListDisSL]
%     - all files are now closed when session is aborted [abortExp]
%     - the photodiode pads are now shown on all the angles of the screen
%     - removed "subject" from filename [outputFile]
%
%

% function BBT_SL_v9

close all;      %close all open figures
clear;          %clear all workspace variables
sca;            %close all screens
% clc;          %clear command window
rng('shuffle')  %reset the random seed

%% Session info

exp_version = 'color'; %set the version of the experiment: 'bw' or 'color'
debug = 1;

if debug == 1 %|| strcmp(getenv('ComputerName'), 'COLLES-D152687')
    trigger = 0;
    projector = 0;
    el.eyelink = 0;
else
    trigger = 1;
    projector = 1;
    el.eyelink = 0;
end

[subject_nr, subject_id, exp_session, exp_part, practice, screen_cm, view_distance] = sessionInfo(debug);
exp_dir = cd;

%% Screen init 1

% Colors
black = 0;
white = .8;

background_color = .102;
pad_color = white;
text_color = white;

% Screen init 1
Screen('Preference', 'SkipSyncTests', debug); %set to 0 during the experiment, 1 during testing
AssertOpenGL;
PsychDefaultSetup(2); %default settings, and unit color range
if strcmp(getenv('ComputerName'), 'COLLES-D152687')
    screenid = 1;
else
    screenid = max(Screen('Screens')); %choose the screen with maximum id (secondary display)
end
[win, rect] = PsychImaging('OpenWindow', screenid, background_color);
[scr_width, scr_height] = Screen('WindowSize', win);

%% Eyelink init

if el.eyelink
    el.edf_file = [subject_nr subject_id '_' num2str(exp_session) num2str(exp_part) '.edf']; %edf filename (max 8 chars)
    el.eyelink_key = KbName('E');  %key used to toggle eyelink validation/calculation on or off during experiment.
    el.continue_key = KbName('C'); %skip eye check part in start/end box, inorder to continue exp
    [el.eye_dir] = eyelinkSetup(debug, el.eyelink, exp_dir, el.edf_file);
end

%% Screen init 2

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');

max_priority = MaxPriority(win); 
if ~debug
    HideCursor; %hide the mouse cursor
    Priority(max_priority); %give max priority to psychtoolbox
end

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Flip', win);

ifi = Screen('GetFlipInterval', win); %get the flip interval of the screen (duration of a frame = 1/refresh rate)
frate = Screen('NominalFrameRate', win); %frame rate
% Screen('TextSize', window, 25);

%% Trigger init

[trig_port_address, trig_iohandle, ~] = triggerSetup(trigger);

% trigger bits
STI001 = 1;     STI002 = 2;     STI003 = 4;     STI004 = 8;
STI005 = 16;	STI006 = 32;    STI007 = 64;    STI008 = 128;

% trigger codes
trig_off = 0;
trig_trial = STI001;
trig_ftag = STI002;
trig_search = STI003;

%% Keyboard init

nataboxSetting(projector)

%% Experimental variables

% Display params
degrees_from_fix = 5; % set the distance from the fixation point (i.e., degrees of the visual field) used for stimulus presentation
stim_size_deg = 6; % stimulus size in degrees
[ecc_xy, degToPixels, stim_size] = degreeToPixels(rect, screen_cm, view_distance, degrees_from_fix, stim_size_deg);
fix_size = degToPixels * .33; %size of the fixation dot
pad_size = 25; %size of the corner pads

% Design params
if debug == 1
    block_nr = 1;
else
    block_nr = 30; %2 times 30 blocks, for a total of 60 blocks per session
end

if practice == 1
    [condition_list, trial_nr] = trialListPractice;
else
    [condition_list, trial_nr] = trialListDisSL(block_nr, subject_nr, exp_session);
end

tlocs_list = condition_list(:,1); %target location (1:upper-left, 2:lower-left, 3:upper-right, 4:lower-right)
dlocs_list = condition_list(:,2); %distractor location (0:absent, 1-4:as above)
toris_list = condition_list(:,3); %target orientation (1:leftward, 2:rightward)
block_list = condition_list(:,4); %block number

% Timing params
dur_fix = .5; % fix display duration
dur_placeholder = 1.5; % placeholder display duration
dur_search = .3; % search display duration
dur_responsetimeout = 2; % max search display duration
dur_itiM = 1; %mean of the inter-trial interval
dur_itiD = .5; %deviation of the inter-trial interval
dur_iti_list = dur_itiM + (2 * rand(trial_nr, 1) - 1) * dur_itiD; dur_iti_list = round(dur_iti_list, 2);
dur_bbt = dur_placeholder + dur_search;

% Stimulus params (gabor patches)
phases = linspace(0,360,5); phases = phases(1:end-1); %phase of the gabors sine grating in degrees
freq_degree = 2; freq = freq_degree/degToPixels; %%spatial frequency in cycle/degree and then in cycles/pixel
sc = stim_size / 6.66; %spatial constant (sigma) of the gaussian hull function of the gabor (i.e., standard deviation)

[gaborID, gaborRect] = makeGabor(win, exp_version, stim_size, freq, sc, phases);

tilt_tar = [345 15]; %orientation(s) (in degrees) of the target
tilt_dis = 90;

tilt_tar_list = tilt_tar(toris_list); %assign a random orientation to the target for each trial
tilts_list = zeros(trial_nr,4); %make a matrix of orientation trial*stimulus location
for i = 1:trial_nr
    tilts_list(i,tlocs_list(i)) = tilt_tar(toris_list(i));
    if dlocs_list(i) %if d_loc is not zero (i.e., it is present)...
        tilts_list(i,dlocs_list(i)) = tilt_dis; %...assign tilt_dis to d_loc
    end
end

if strcmp(exp_version, 'color') %assign a random color to the target for each trial
    col_tar = ['g', 'r'];
    tcol_list = repmat(1:2, trial_nr/2, 1);
    tcol_list = tcol_list(randperm(numel(tcol_list)));
end

% Location params
[xy, rect_ul, rect_ur, rect_ll, rect_lr] = propixxXY(gaborRect, scr_width, scr_height, ecc_xy);

% Broadaband tagging params
if strcmp(exp_version, 'bw')
    frate_incr = 12; %increment in frate generated by the Propixxx projector (4 quadrants * 3 color channels = 12)
elseif strcmp(exp_version, 'color')
    frate_incr = 4; %increment in frate generated by the Propixxx projector (4 quadrants = 4)
end
frate_rapidmode = frate * frate_incr;

% [tag_list, dur_bbt_propixx] = bbTaggingSignal(exp_version, trial_nr, frate_rapidmode, dur_bbt); %generate signals
load tag_list;
dur_bbt_propixx = dur_bbt * frate_rapidmode;

dur_fix_propixx = dur_fix * frate;
dur_bbt_propixx = dur_bbt_propixx / frate_incr;
dur_placeholder_propixx = dur_bbt_propixx * dur_placeholder/dur_bbt;
dur_search_propixx = dur_bbt_propixx * dur_search/dur_bbt;

[pad_frame, rect_DLpad, rect_DRpad, rect_ULpad, rect_URpad] = cornerPads(xy, pad_size, pad_color); %draw photodiode frames
padID = Screen('MakeTexture', win, pad_frame);

% Response params
resp_key = zeros(trial_nr, 1);
resp_cor = zeros(trial_nr, 1);
resp_rt = zeros(trial_nr, 1);

if projector == 1
  nata_keys = [4 7];
  toris_list = nata_keys(toris_list);
end

% Get date and time
datetime = datestr(now);
    
% Log file
var_list = 'subject_nr\t subject_id\t datetime\t session_count\t part_count\t block_count\t trial_count\t target_loc\t distractor_loc\t target_col\t target_ori\t target_ori_key\t correct\t response\t RT'; % list of the variables present in output header
[outputfile] = outputFile(subject_nr, subject_id, exp_session, exp_part, var_list);

%% Open projector

% Get projector mode
if  strcmp(exp_version, 'bw')
    proj_mode = 5; %b/w 1440 Hz
elseif strcmp(exp_version, 'color')
    proj_mode = 2; %color 480 Hz
end

% Open projector
propixxSetup('open', projector, proj_mode);

%% Eye link calibration & validation

% Reset projector
propixxSetup('reset', projector, proj_mode);

% Set parameters, start and calibrate eyelink
if el.eyelink
    [el.defaults] = eyelinkStart(el.edf_file, screenid, rect);
end

%% Trial run

% Set projector
propixxSetup('set', projector, proj_mode);

% Start eyelink
if el.eyelink
    Eyelink('Message', 'Experiment start'); %experiment start message to eyelink
    Eyelink('Message', ['Block ' int2str(exp_part) ' start']);
end

% Draw initial message
propixxDrawText(win, 'Press any key to begin', xy, text_color)
if debug == debug
  Screen('DrawTextures', win, padID, [], rect_DLpad, [], [], [], []);
  Screen('DrawTextures', win, padID, [], rect_DRpad, [], [], [], []);
  Screen('DrawTextures', win, padID, [], rect_ULpad, [], [], [], []);
  Screen('DrawTextures', win, padID, [], rect_URpad, [], [], [], []);
end
Screen('Flip', win);
KbWait;
Screen('Flip', win);
WaitSecs(1);

for triali = 1:trial_nr
   
    % Set variable values
    t_locL = tlocs_list(triali);
    if strcmp(exp_version, 'color')
        tcolL = tcol_list(triali);
    end
    d_locL = dlocs_list(triali);
    t_tiltL = tilt_tar_list(triali);
    
    % Randomize gabor phases
    phasesL = randperm(length(phases)); 
    
    % Select gabor textures
    gaborID_p = zeros(4,1);
    gaborID_s = zeros(4,1);
    if strcmp(exp_version, 'color')
        if tcolL == 1 %green
            gaborID_p(1) = gaborID(phasesL(1)); %placeholder array
            gaborID_p(2) = gaborID(phasesL(2));
            gaborID_p(3) = gaborID(phasesL(3));
            gaborID_p(4) = gaborID(phasesL(4));
            
            gaborID_s(1) = gaborID(phasesL(1)); %search array
            gaborID_s(2) = gaborID(phasesL(2));
            gaborID_s(3) = gaborID(phasesL(3));
            gaborID_s(4) = gaborID(phasesL(4));
            if d_locL
                gaborID_s(d_locL) = gaborID(phasesL(d_locL) + 4);
            end
        else
            gaborID_p(1) = gaborID(phasesL(1) + 4); %placeholder array
            gaborID_p(2) = gaborID(phasesL(2) + 4);
            gaborID_p(3) = gaborID(phasesL(3) + 4);
            gaborID_p(4) = gaborID(phasesL(4) + 4);
            
            gaborID_s(1) = gaborID(phasesL(1) + 4); %search array
            gaborID_s(2) = gaborID(phasesL(2) + 4);
            gaborID_s(3) = gaborID(phasesL(3) + 4);
            gaborID_s(4) = gaborID(phasesL(4) + 4);
            if d_locL
                gaborID_s(d_locL) = gaborID(phasesL(d_locL));
            end
        end
        tcolL = col_tar(tcolL);
    elseif strcmp(exp_version, 'bw')
        gaborID_p(1) = gaborID(phasesL(1)); %placeholder array
        gaborID_p(2) = gaborID(phasesL(2));
        gaborID_p(3) = gaborID(phasesL(3));
        gaborID_p(4) = gaborID(phasesL(4));

        gaborID_s(1) = gaborID(phasesL(1)); %search array
        gaborID_s(2) = gaborID(phasesL(2));
        gaborID_s(3) = gaborID(phasesL(3));
        gaborID_s(4) = gaborID(phasesL(4));
        
        tcolL = exp_version;
    end
    
    % Select frequency tagging signal
    tagL = cell2mat(tag_list(1,triali));
    tagR = cell2mat(tag_list(2,triali));
    
    % Set triggers
    d_locTemp = 1:4; d_locTemp(d_locTemp == t_locL) = 0;   %set distractor absent equal to target location
    d_locTemp = find(d_locTemp == d_locL);
    C1 = ~rem(d_locTemp, 2);    %distractor location: top/bottom (0:up, 1:down)
    C2 = d_locTemp > 2;         %distractor location: side (0:left, 1:right)
    C3 = ~rem(t_locL, 2);       %target location: top/bottom (0:up, 1:down)
    C4 = t_locL > 2;            %target location: side (0:left, 1:right)
    C5 = t_tiltL == 15;         %target orientation (0:left[345], 1:right[15])
    trig_condi = C1 * STI004 + C2 * STI005 + C3 * STI006 + C4 * STI007 + C5 * STI008;
    
    % Send initial message to eyelink
    if el.eyelink %send trial start trigger to eyelink
        Eyelink('Message', [int2str(exp_part) ' ' int2str(triali) ' Distractor location: ' int2str(d_locL) ' Target location: ' int2str(t_locL) ' Target orientation: ' int2str(t_tiltL)]);
    end
    
    % Send start trial trigger
    triggerSend(trigger, trig_iohandle, trig_port_address, trig_trial);
    WaitSecs(.05);
    triggerSend(trigger, trig_iohandle, trig_port_address, trig_off);
    
    % Draw fixation point
    vbl = GetSecs;
    for di = 1:dur_fix_propixx
        Screen('DrawDots', win, xy, fix_size, text_color, [], 1);
        Screen('DrawDots', win, xy, fix_size/4, background_color, [], 1);
        if di == 1
            [trig_time, trig_sent] = triggerSend(trigger, trig_iohandle, trig_port_address, trig_trial);
        end
        vbl = Screen('Flip', win, vbl + .5 * ifi);
        if trig_sent && GetSecs > (trig_time + .05)
            [~, trig_sent] = triggerSend(trigger, trig_iohandle, trig_port_address, trig_off);
        end
    end
    
    % Draw the placeholder and the search array
    for di = 1:dur_bbt_propixx
        tagL_frame =  tagL(:, di);
        tagR_frame =  tagR(:, di);
        if di <= dur_placeholder_propixx %placeholder array
          for i = 1:4
            Screen('DrawTexture', win, gaborID_p(1), [], rect_ul(:,i), [], [], [], .5, [], 0);
            Screen('DrawTexture', win, gaborID_p(2), [], rect_ll(:,i), [], [], [], tagL_frame(i,:), [], 0);
            Screen('DrawTexture', win, gaborID_p(3), [], rect_ur(:,i), [], [], [], .5, [], 0);
            Screen('DrawTexture', win, gaborID_p(4), [], rect_lr(:,i), [], [], [], tagR_frame(i,:), [], 0);
            Screen('DrawDots', win, xy(:,i), fix_size, text_color, [], 1);
            Screen('DrawDots', win, xy(:,i), fix_size/4, background_color, [], 1);
            Screen('DrawTexture', win, padID, [], rect_DLpad(:,i), [], [], [], tagL_frame(i,:));
            Screen('DrawTexture', win, padID, [], rect_DRpad(:,i), [], [], [], tagR_frame(i,:));
            Screen('DrawTexture', win, padID, [], rect_ULpad(:,i), [], [], [], tagL_frame(i,:));
            Screen('DrawTexture', win, padID, [], rect_URpad(:,i), [], [], [], tagR_frame(i,:));
          end
        else
          for i = 1:4 %search array
            Screen('DrawTexture', win, gaborID_s(1), [], rect_ul(:,i), tilts_list(triali,1), [], [], .5, [], 0);
            Screen('DrawTexture', win, gaborID_s(2), [], rect_ll(:,i), tilts_list(triali,2), [], [], tagL_frame(i,:), [], 0);
            Screen('DrawTexture', win, gaborID_s(3), [], rect_ur(:,i), tilts_list(triali,3), [], [], .5, [], 0);
            Screen('DrawTexture', win, gaborID_s(4), [], rect_lr(:,i), tilts_list(triali,4), [], [], tagR_frame(i,:), [], 0);
            Screen('DrawDots', win, xy(:,i), fix_size, text_color, [], 1);
            Screen('DrawDots', win, xy(:,i), fix_size/4, background_color, [], 1);
            Screen('DrawTexture', win, padID, [], rect_DLpad(:,i), [], [], [], tagL_frame(i,:));
            Screen('DrawTexture', win, padID, [], rect_DRpad(:,i), [], [], [], tagR_frame(i,:));
            Screen('DrawTexture', win, padID, [], rect_ULpad(:,i), [], [], [], tagL_frame(i,:));
            Screen('DrawTexture', win, padID, [], rect_URpad(:,i), [], [], [], tagR_frame(i,:));
          end
        end
        if di == 1
            [trig_time, trig_sent] = triggerSend(trigger, trig_iohandle, trig_port_address, (trig_condi + trig_ftag));
        elseif di == dur_placeholder_propixx + 1
            [trig_time, trig_sent] = triggerSend(trigger, trig_iohandle, trig_port_address, trig_search);
        end
        vbl = Screen('Flip', win, vbl + .5 * ifi);
        if trig_sent && GetSecs > (trig_time + .05)
            [~, trig_sent] = triggerSend(trigger, trig_iohandle, trig_port_address, trig_off);
        end
        if di == dur_placeholder_propixx + 1
          start_time = vbl;
        end
    end
    
    % Draw blank screen
    Screen('DrawDots', win, xy, fix_size, text_color, [], 1);
    Screen('DrawDots', win, xy, fix_size/4, background_color, [], 1);
    Screen('Flip', win);
    
    % Get keyboard response
    [resp_key(triali), resp_rt(triali)] = nataboxResponse(dur_responsetimeout, start_time, debug, el);
    resp_keyL = resp_key(triali);
    resp_rtL = resp_rt(triali);
    
    % ITI and errror message
    cor_resp = toris_list(triali); %get correct response
    resp_cor(triali) = resp_key(triali) == cor_resp;
    resp_corL = resp_cor(triali);
    
    if ~resp_cor(triali)
        propixxDrawText(win, 'ERROR', xy, text_color)
    end
    dur_iti = dur_iti_list(triali);
    Screen('Flip', win);
    WaitSecs(dur_iti);

    % Break
    if triali ~= trial_nr && rem(triali, 12*10) == 0
        propixxDrawText(win, 'Take a break. Press any key when you are ready to continue', xy, text_color)
        Screen('Flip', win);
        KbWait;
    end

    % Save results to output file
    blockL = block_list(triali); %get current block number
    fprintf(outputfile, '\n%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d',...
        subject_nr, subject_id, datetime, exp_session, exp_part, blockL, triali, ...
        t_locL, d_locL, tcolL, t_tiltL, cor_resp, resp_corL, resp_keyL, resp_rtL);

end

% End message
propixxDrawText(win, 'The experiment is finished. Thanks for your participation', xy, text_color)
Screen('Flip', win);
KbWait;

%% Closing

% reset projector
propixxSetup('reset', projector, proj_mode);
propixxSetup('close', projector, proj_mode);

% stop eyelink
if el.eyelink
    Eyelink('Message', 'END OF SESSION');
    eyelinkStop(el);
end

save(sprintf('Data/%s_%s_s%sp%s', num2str(subject_nr), subject_id, num2str(exp_session), num2str(exp_part))); %save workspace variables to file
Screen('Close', win); %closes the window
fclose('all'); %closes all files (e.g., the log file)
Priority(0); %removes the priority
ShowCursor; %shows the cursor
KbQueueRelease() %clean up the KbQueue

% end

%% FUNCTIONS

% sessionInfo function
function [subject_nr, subject_id, exp_session, exp_part, practice, screen_cm, view_distance] = sessionInfo(debug)
% Show a window and ask session information

    prompt = {'Participant Nr:', 'Participant ID:', 'Session', 'Part', 'Practice', 'Screen size: ', 'Viewing distance:'}; %description of fields
    if debug || strcmp(getenv('ComputerName'), 'COLLES-D152687')
        defaults = {'999','none','0','0','0','50.8','57'}; %default values for debuging
    else
        defaults = {'','','','','0','70','145'};%default values 
    end
    answer = inputdlg(prompt, 'Participant', 1.2, defaults); %opens dialog

    % Check if there is already a subject with the insert id
    if isempty(answer)
        error('Program aborted.')
    elseif isempty(answer{1,:}) || isempty(answer{2,:}) || isempty(answer{3,:})
        error('Program aborted. Session info required.')
    elseif exist(['Data/subject_' num2str(answer{1,:}) '_' answer{2,:} '_part' num2str(answer{3,:}) '.txt'], 'file') == 2
        fclose('all'); % close all files
        button = questdlg('Overwrite subject data?', '','Yes','No','No'); % ask to overwrite subject data
        if strcmp(button, 'Yes')
            delete(['Data/subject_' num2str(answer{1,:}) '_' answer{2,:} '_part' num2str(answer{3,:}) '.txt']); % if yes, delete old data file and continue the experiment
        else
            error('Program aborted. Participant ID already used.') % if no, abort the experiemnt and show an error message
        end
    end

    % Assign values to the corrisponding variables
    subject_nr = str2double(answer{1,:}); %subject number
    subject_id = answer{2,:}; %subject ID (from participant PC)
    exp_session = str2double(answer{3,:}); %experiment session (1 or 2)
    exp_part = str2double(answer{4,:}); %experiment part (1 or 2)
    practice  = str2double(answer{5,:}); %practice (set to 1 to run practice)
    screen_cm = str2double(answer{6,:}); %screen size in cm
    view_distance = str2double(answer{7,:}); %viewing distance
    
end

% outputFile function
function [outputfile] = outputFile(subject_nr, subject_id, exp_session, exp_part, varList)
% Create output file

    % Check and make data directory
    if ~exist('data', 'dir') % isfolder('data')
        mkdir Data
    end
    
    % Set data folder
    dataFolder = 'Data'; 
    
    % Create the output file
    outputfile = fopen([dataFolder '/' num2str(subject_nr) '_' subject_id '_s' num2str(exp_session) 'p' num2str(exp_part) '.txt'],'a');
    fprintf(outputfile, varList); %print variable list

end

% degreesToPixels function
function [eccXY, degToPixels, stimSizeXY] = degreeToPixels(rect, screen_cm, viewDistance, degreesFromFix, stimSizeDeg)
% Translate the size of one degree of visual filed in pixels.
% 'viewDistance' is the distance between the participant's eyes (mesueres
% from the nasion) and the center of the screen;
% 'degreesFromFix' is the stimulus eccentricity value (i.e., the distance
% between each stimulus and the central fixation point/cross) ion degrees.
    
    degToCms = 2*viewDistance*tan(pi/360); %get the measure in cm of one degree based on the viewing distance
    
    if 1 %implement this conditional if you want to use a normal mode presentation
        cmToPixels = screen_cm/(rect(3)/2); %get the measure of one pixel in cm (cm/px)
    %else
    %    cmToPixels = screen_cm/rect(3); %get the measure of one pixel in cm (cm/px)
    end
    
    degToPixels = round(degToCms/cmToPixels); %get the number of pixels in one degree
    ecc_degToPixels = degreesFromFix * degToPixels; %transform eccentricity from degrees to pixels 
    eccXY = round(ecc_degToPixels/sqrt(2)); %get the cathetus (i.e., the x - and y - value) from the hypotenuse (i.e., the eccentricity)
    %Note the here there are four location. Thus, the x and y values are
    %identical. With more locations (e.g, 6), the two values must be computed
    %separately.
    
    stimSizeXY = round(stimSizeDeg * degToPixels); %get stimulus size in pixels.
    %Again, x and y have the same value. If the stimuli are not equilateral, 
    %compute x and y separately.
    
end

% trialListPractice function
function [TrialL, TrialNr] = trialListPractice
% Build practice trial table
    
    % Trial table
    TLOC = repmat(1:4, [5 1]);
    TLOC = TLOC(:);
    DLOC = repmat(0:4, [1 4]);
    DLOC = DLOC(:);
    DESIGN = [TLOC DLOC];
    
    % Remove all trials with TLOC equal to DLOC
    DESIGN = DESIGN(~bsxfun(@eq, DESIGN(:,1), DESIGN(:,2)),:);

    % Randomize trials
    blocks = 1;
    TrialL = zeros(size(DESIGN,1)*blocks,size(DESIGN,2));
    for bi = 1:blocks
        DESIGN_RAND = Shuffle(DESIGN,2);
        TrialL((1+(size(DESIGN,1)*(bi-1))):(size(DESIGN,1)*bi),:) = DESIGN_RAND;
    end
    
    TrialNr = size(TrialL, 1);
    
    % Set target orientations randomly (not counterbalanced)
    toris = repmat(1:2, TrialNr/2, 1);
    TrialL(:,3) = toris(randperm(numel(toris)));
    
    % Set block number (0 = practice)
    TrialL(:,4) = 0;
    
end

% trialListForSL function
function [TrialL, TrialNr] = trialListDisSL(block_number, subject_nr, exp_session)
% Build trial table with unbalanced distractor presentiation across
% locations.
    
    % Set target and distractor params
    %block_number = 4; %number of blocks (for debugging)
    tloc_disabs = [1 2 3 4]; %target location in distractor absent trials
    tloc_dispre = [1 1 2 2 3 3 4 4]; %target location in distractor present trials
    dloc_disabs = [0 0 0 0]; %distractor location in distractor absent trials
    
    % Select probability-hemifield association based on subject_nr and exp_session
    if mod(subject_nr, 2) == 1
        if exp_session < 2
            dloc_dispre = [1 1 1 2 2 2 3 4]; %distractor location in distractor present trials
        else
            dloc_dispre = [1 2 3 3 3 4 4 4];
        end
    else
        if exp_session < 2
            dloc_dispre = [1 2 3 3 3 4 4 4];
        else
            dloc_dispre = [1 1 1 2 2 2 3 4];
        end
    end
    
    % Generates random pairs of target and distractor locations
    TrialL = zeros(length(tloc_disabs)*2*block_number,4);
    
    for blocki = 1:block_number
        tloc_dispre_temp = tloc_dispre;
        
        while any(bsxfun(@eq, tloc_dispre_temp, dloc_dispre)) == 1
            tloc_dispre_temp = tloc_dispre(randperm(length(tloc_dispre)));
        end
        
        TLOC = [tloc_disabs tloc_dispre_temp];
        TLOC = TLOC(:);
        DLOC = [dloc_disabs dloc_dispre];
        DLOC = DLOC(:);
        
        DESIGN = [TLOC DLOC];
        DESIGN_RAND = Shuffle(DESIGN,2);
        
        TrialL((1+(size(DESIGN,1)*(blocki-1))):(size(DESIGN,1)*blocki),1:size(DESIGN,2)) = DESIGN_RAND;
        TrialL((1+(size(DESIGN,1)*(blocki-1))):(size(DESIGN,1)*blocki),4) = blocki;
    end
    
    TrialNr = size(TrialL, 1);
    
    % Set target orientations randomly (not counterbalanced)
    toris = repmat(1:2, TrialNr/2, 1);
    TrialL(:,3) = toris(randperm(numel(toris)));

end

% triggerSetup function
function [port_address, iohandle, status] = triggerSetup(trigger)

    if trigger
        port_address = hex2dec('BFF8'); %check this port address!
        iohandle = io64;
        status = io64(iohandle); % initialize the interface to the inpoutx64 system driver
        io64(iohandle, port_address, 0); %send 0 trigger (reset all pins)
    else
        port_address = [];
        iohandle = [];
        status = [];
    end
    
end

% triggetSend function
function [time, sent]= triggerSend(trigger, iohandle, port_address, trig_code)

    if trigger
        io64(iohandle, port_address, trig_code); %send trigger to MEG
    end

    time = GetSecs;
    if trig_code
        sent = 1;
%         if ~trigger
%             display(trig_code)
%         end
    else
        sent = 0;
    end
    
end

% nataboxSetting function
function nataboxSetting(projector)

    device = 0;
    KbName('UnifyKeyNames'); % for easy use of Keyboard keys
    
    if projector == 0
        response_keys = [KbName('1') KbName('2')]; % debugging allowed keys
    else
        response_keys = [KbName('4$') KbName('7&')]; %These are the left and right index fingers of the (5-button) NATA boxes
    end
    
    active = [response_keys, KbName('END')]; %by default, the end key is allowed
    keylist=zeros(1,256); %Set all keys to zero (ignore)
    keylist(active)=1; %set active keys to 1 (listen)
    KbQueueCreate(device,keylist);%%Create queue, this is a time consuming operation (relatively), do while non-time critical
%     KbQueueStart(); %Start listening
    
end

% nataboxResponse fnction
function [response_key, response_time] = nataboxResponse(response_timeout, start_time, debug, el)

    KbQueueStart(); %Start listening
    
    if nargin < 2
        start_time = GetSecs; %to get reaction times relative to some event (e.g. stimulus onset)
    end

    KbQueueFlush();% clear all keyboard presses so far. Basically, you want the responses after this point 

    %listen to reponses for a fixed amount of time, using KbQueueCheck (similar to KbCheck)
    while (GetSecs - start_time) < response_timeout
        
        response_key = 0;
        response_time = 0;
        
        [pressed, firstpress] = KbQueueCheck(); %check response, return whether pressed, and first press timestamp
        
        if any(pressed) %exit loop once a response is recorded
            keyCode = find(firstpress);
            if keyCode == KbName('END') %end key quits the experiment
              abortExp(debug, el)
            else
              keyName = KbName(keyCode);
              response_key = str2double(keyName(1));
              rTime = firstpress(keyCode);
              response_time = (rTime-start_time) * 1000;
              break
            end
        end
        
    end
    
end

% propixxSetup function
function propixxSetup(command, projector, mode)
% Open, set, reset and close Propixx projector

    if projector == 1
      if strcmp(command, 'open')
        Datapixx('Open');
      elseif strcmp(command, 'set')
        Datapixx('SetPropixxDlpSequenceProgram', mode); % 2= color 480 Hz, 5= b/w 1440 Hz
        Datapixx('RegWrRd');
      elseif strcmp(command, 'reset')
        Datapixx('SetPropixxDlpSequenceProgram', 0); % default
        Datapixx('RegWrRd');
      elseif strcmp(command, 'close')
        Datapixx('Close');
      else
        fprintf(1, 'Propixx command ''%s'' is not defined.\n', command);
        return
      end 
    end
  
end

% propixxXY function
function [xy, rect_ul, rect_ur, rect_ll, rect_lr] = propixxXY(stim_rect, scr_width, scr_height, ecc_xy)
% Divide the screen in four quadrants and return the coordinates of the
% center of each quadrant (xy) and the rects of the four location used for
% stimulus presentation (upper-let, upper-right, lower-left, lopwer-right)

    % centers
    dstrect_1 = CenterRectOnPoint(stim_rect, 1 * scr_width / 4, 1 * scr_height / 4);
    dstrect_2 = CenterRectOnPoint(stim_rect, 3 * scr_width / 4, 1 * scr_height / 4);
    dstrect_3 = CenterRectOnPoint(stim_rect, 1 * scr_width / 4, 3 * scr_height / 4);
    dstrect_4 = CenterRectOnPoint(stim_rect, 3 * scr_width / 4, 3 * scr_height / 4);
    dstrect = {dstrect_1, dstrect_2, dstrect_3, dstrect_4};

    % center, left and right below corner
    xy = zeros(2, 4); 
    rect_ul = zeros(4, 4); 
    rect_ur = zeros(4, 4);
    rect_ll = zeros(4, 4); 
    rect_lr = zeros(4, 4);
    for i = 1:4
      [x, y] = RectCenter(dstrect{i});
      xy(:, i) = [x, y];
      ul = CenterRectOnPoint(dstrect{i}, x - ecc_xy, y - ecc_xy);
      rect_ul(:, i) = ul;
      ur = CenterRectOnPoint(dstrect{i}, x + ecc_xy, y - ecc_xy);
      rect_ur(:, i) = ur;
      ll = CenterRectOnPoint(dstrect{i}, x - ecc_xy, y + ecc_xy);
      rect_ll(:, i) = ll;
      lr = CenterRectOnPoint(dstrect{i}, x + ecc_xy, y + ecc_xy);
      rect_lr(:, i) = lr;
    end

end

% propixxDrawText function
function propixxDrawText(win, text, xy, text_color)

    for i = 1:4
        x = xy(1,i) - 12 * length(text) / 2;
        y = xy(2,i);
%         Screen('DrawText', win, text, x, y, text_color);
        DrawFormattedText(win, text, x, y, text_color);
    end
    
end

% makeGabor function
function [gaborID, gaborRect] = makeGabor(win, exp_version, stim_size, freq, sc, phases)

%     %debuging
%     stim_size = 152;
%     sc = 25;
%     freq = .1;
%     phases = linspace(0,360,5);
%     phases = phases(1:end-1);
%     exp_version = 'bw'; %'color'
    
    stim_size_half = stim_size/2;
    
    % colors
    Color_1 = [0.33, 0.72, 0.33]; %green
    Color_2 = [0.87, 0.33, 0.33]; %red
    Color = [Color_1; Color_2];
    
    % phases
    phi = zeros(1, numel(phases));
    for i = 1:numel(phases)
        phi(i) = (phases(i) / 180) * pi;
    end

    % make sinusoid
    [gX, gY] = meshgrid(1:stim_size, 1:stim_size);

    Xs = zeros(stim_size, stim_size, numel(phases));
    for i = 1:numel(phases)
        Xs(:,:,i) = sin(2 * pi * freq * (gX - stim_size_half) + phi(i));
    end
    Xs = (Xs + 1) / 2; %normalization
    
    % make gaussian
    G = exp(-((gX - stim_size_half) .^ 2 + (gY - stim_size_half) .^ 2) ./ (2 * sc ^ 2));
    
    % make circular mask
    mask = (gX - stim_size_half) .^ 2 + (gY - stim_size_half) .^ 2 <= (stim_size_half) ^ 2;
    
    % combine gaussian and mask
    Gmask = G .* mask;
    
    % make gabor
    if strcmp(exp_version, 'bw')
        gabor = cell(numel(phases),1);
        for i = 1:numel(phases)
            Y = zeros(stim_size, stim_size, 4);
            Y(:, :, 1) = Xs(:,:,i);	% red
            Y(:, :, 2) = Xs(:,:,i);	% green
            Y(:, :, 3) = Xs(:,:,i);	% blue
            Y(:, :, 4) = Gmask;     % alpha
            gabor{i} = uint8(Y * 255);
        end
    elseif strcmp(exp_version, 'color')
        gabor = cell(numel(phases),2);
        for i = 1:numel(phases)
            for j = 1:2
                Y = zeros(stim_size, stim_size, 4);
                Y(:, :, 1) = Xs(:,:,i) *(Color(j,1));	% red
                Y(:, :, 2) = Xs(:,:,i) *(Color(j,2));	% green
                Y(:, :, 3) = Xs(:,:,i) *(Color(j,3));	% blue
                Y(:, :, 4) = Gmask;                     % alpha
                gabor{i, j} = uint8(Y * 255);
            end
        end
    end
        
    % make gabor textures and get ids
    if strcmp(exp_version, 'bw')
        gaborID = zeros(1, numel(phases));
        for i = 1:numel(phases)
          gaborID(i) = Screen('MakeTexture', win, gabor{i});
        end
    elseif strcmp(exp_version, 'color')
        gaborID = zeros(numel(phases),2);
        for i = 1:numel(phases)
          for j = 1:size(Color,1)
            gaborID(i, j) = Screen('MakeTexture', win, gabor{i, j});
          end
        end
        gaborID = gaborID(:);
    end
    
    % get rects
    gaborRect = Screen('Rect', gaborID(1));
    
end

% eyelinkSetup
function [eye_dir] = eyelinkSetup(debug, eyelink, exp_dir, edf_file)
    
    override = 0;

    if eyelink
        
        %add eyelink script folder (should be in main experiment folder)
        addpath([exp_dir filesep 'Eyelink']);
        
        %make directory if it doesn't already exist (local computer)
        eye_dir = [exp_dir filesep 'Eyelink' filesep ];
        if ~exist(eye_dir, 'dir')
            mkdir(eye_dir);
        end
        
        %check whether files already exist for this subject/session
        if exist([exp_dir filesep 'Eyelink' filesep 'Data' filesep  edf_file '.edf'],'file')>0
            cont = input('Warning! Eyelink file will be overwritten, do you want to continue? (y/n) ','s');
            if cont == 'n'
                error('Session aborted')
            end
        end
        
%         rect = [0 0 resx resy]; %% needed for the el_Set_Params function
        
    else % when eyelink is off
        if ~debug %is this is real experiment time, eyelink should be on
            if override
                warning('Eyelink not in use, continuing anyway...') %#ok<UNRCH>
            else
                error('Eyelink not selected!')
            end
        end
    end

end

% eyelinkParams function
function [defaults] = eyelinkParams(defaults, rect)
%Custom parameters for eyelink

    el = defaults;

    el.eye_used                = 'LEFT_EYE';
    el.calibrationtargetsize   = 1;
    el.calibrationtargetwidth  = 0.5;
    el.targetbeep              = 0;
    el.feedbackbeep            = 0;
    el.displayCalResults       = 1;
    el.eyeimagesize            = 50;  % percentage of screen
    %el.backgroundcolour        = 127.5; GrayIndex(el.window);

    disp('Updating Parameters')
    EyelinkUpdateDefaults(el);

    defaults = el;

    % make sure we're still connected.
    if Eyelink('IsConnected')~=1
        warning('eyelink is not connected! restart the tracker');
        Eyelink('Shutdown'); %shutdown Eyelink:
        el.online = 0;
        return;
    end

    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT'); 

    % This Command is crucial to map the gaze positions from the tracker to
    % screen pixel positions to determine fixation
    Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld',  rect(1)-rect(1), rect(2)-rect(2), rect(3)-rect(1), rect(4)-rect(2));
    Eyelink('message','DISPLAY_COORDS %ld %ld %ld %ld',         rect(1)-rect(1), rect(2)-rect(2), rect(3)-rect(1), rect(4)-rect(2));

    % Use Psychophysical setting
    Eyelink('Command', 'recording_parse_type = GAZE');
    Eyelink('Command', 'saccade_velocity_threshold = 22');
    Eyelink('Command', 'saccade_acceleration_threshold = 3800');
    Eyelink('Command', 'saccade_motion_threshold = 0.0');
    Eyelink('Command', 'saccade_pursuit_fixup = 60');
    Eyelink('Command', 'fixation_update_interval = 0');

    % Other tracker configurations

    % these might crash:
    Eyelink('Command', 'heuristic_filter = 0');
    Eyelink('Command', 'pupil_size_diameter = YES');

    % use 9 point calibration (Default)
    Eyelink('Command', 'calibration_type = HV9');
    %Eyelink('Command', 'calibration_type = HV13');

    Eyelink('Command', 'generate_default_targets = YES');
    Eyelink('Command', 'enable_automatic_calibration = YES');
    Eyelink('Command', 'automatic_calibration_pacing = 1000');
    Eyelink('Command', 'binocular_enabled = NO');
    Eyelink('Command', 'use_ellipse_fitter = NO');
    Eyelink('Command', 'sample_rate = 1000');
    %Eyelink('Command', 'elcl_tt_power = %d', 2); % illumination, 1 = 100%, 2 = 75%, 3 = 50%

    switch el.eye_used
        case 'RIGHT_EYE'
            Eyelink('Command', 'file_event_filter = RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,INPUT');
            Eyelink('Command', 'link_event_filter = RIGHT,FIXATION,FIXUPDATE,SACCADE,BLINK,MESSAGE,INPUT');
        case  'LEFT_EYE'
            Eyelink('Command', 'file_event_filter = LEFT,FIXATION,SACCADE,BLINK,MESSAGE,INPUT');
            Eyelink('Command', 'link_event_filter = LEFT,FIXATION,FIXUPDATE,SACCADE,BLINK,MESSAGE,INPUT');
    end

    Eyelink('Command', 'file_sample_data  = GAZE,GAZERES,HREF,PUPIL,AREA,STATUS,INPUT');

end

% eyelinkStart function
function [defaults] = eyelinkStart(edf_file, screenid, rect)
% Open screen for calibration, calibrate and start recording

    % STEP 1
    % Open a graphics window on the main screen
    % using the PsychToolbox's Screen function.    
    % use the shrunk version of the window
    %window=Screen('OpenWindow', cfg.screenNumber, [] ,cfg.el_rect);
    win = Screen('OpenWindow', screenid);


    % STEP 2
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    % Psychtoolbox defaults function
    defaults = EyelinkInitDefaults(win);

    % Disable key output to Matlab window:
    %ListenChar(2);

    % STEP 3
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit
        fprintf('Eyelink Init aborted.\n');         
        Eyelink('Shutdown');%shutdown Eyelink:
        sca;% close window:
        ListenChar(0);%restore keyboard output to Matlab:
        return;
    else
        disp('Eyelink initizalized')
    end

    % open file to record data to
    disp('Opening EDF file');   
    status = Eyelink('Openfile', edf_file);

    if ~status
        disp('EDF file opened on Eyelink computer')
    else
        error(['Could not open EDF file on Eyelink computer, error: ' int2str(status)])
    end

    % set custom parameters
    disp('Setting parameters')
    [defaults] = eyelinkParams(defaults, rect);

    % STEP 4 
    % Calibrate the eye tracker
    disp('Starting calibration')
    EyelinkDoTrackerSetup(defaults);

    % do a final check of calibration using driftcorrection
    %     EyelinkDoDriftCorrection(el);

    % STEP 5
    % start recording eye position
    disp('Start recording')
    Screen('Close', win);
    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    WaitSecs(0.1);
    % mark zero-plot time in data file
    disp('Sending message')
    Eyelink('Message', 'SYNCTIME');

%     sca
    ListenChar(0);

end

% eyelinkStop function
function eyelinkStop(el)

% STEP 7
% finish up: stop recording eye-movements,
% close graphics window, close data file and shut down tracker
Eyelink('StopRecording');
Eyelink('CloseFile');
% download data file

fprintf('Receiving data file ''%s''\n', el.edf_file);
%status=Eyelink('ReceiveFile');
status=Eyelink('ReceiveFile',el.edf_file,el.eye_dir,1); %transfer file to experiment directory
if status > 0
    fprintf('ReceiveFile status %d\n', status);
end
if 2==exist(el.edf_file, 'file')
    fprintf('Data file ''%s'' can be found in ''%s''\n', el.edf_file, el.eye_dir);
end

% Shutdown Eyelink:
Eyelink('Shutdown');
% Close window:
% sca;
% Restore keyboard output to Matlab:
ListenChar(0);

end

% cornerPads function
function [pad_frame, rect_DLpad, rect_DRpad, rect_ULpad, rect_URpad] = cornerPads(xy, width, pad_color)

    pad_frame = uint8(ones(width, width) * round(255*pad_color));

    rect_DLpad = zeros(4, 4);
    rect_DRpad = zeros(4, 4);
    rect_ULpad = zeros(4, 4);
    rect_URpad = zeros(4, 4);
    for i = 1:4
        y = xy(2,i) + xy(2,2)  - width;
        x = xy(1,i) - xy(1,1);
        rect_DLpad(:, i) = [x, y, x + width, y + width];
        x = xy(1,i) + xy(1,1) - width;
        rect_DRpad(:, i) = [x, y, x + width, y + width];
        
        y = xy(2,i) - xy(2,2);
        x = xy(1,i) - xy(1,1);
        rect_ULpad(:, i) = [x, y, x + width, y + width];
        x = xy(1,i) + xy(1,1) - width;
        rect_URpad(:, i) = [x, y, x + width, y + width];        
    end
    
end

% bbTaggingSignal
function [tag_list, dur_propixx] = bbTaggingSignal(exp_version, trial_nr, frate, dur_tag)

    % %debug
    % dur_tag = 1.8;
    % ifi = 1/120;
    % ifi_dp = ifi / 4; %color version
    % frate = round(1 / ifi_dp);
    % trial_nr = 360;
    
    dur_propixx = dur_tag * frate;
    t = linspace(0, dur_tag, dur_propixx);
    L = length(t);
    %f = linspace(0, frate, dur_propixx);
    
%     % broadband noise
%     fh = 75;
%     fl = 55;
%     [bl, al] = butter(8, fh / (frate / 2), 'low');  
%     [bh, ah] = butter(8, fl / (frate / 2), 'high'); 

    % phase-modulated signal
    [bl, al] = butter(4, 6 / (frate / 2), 'low'); %filter params
    [bh, ah] = butter(4, 2 / (frate / 2), 'high'); 

    tag_list = cell(2,trial_nr);
    %pDelays = -100:1:100;
    %xtag1 = zeros(length(pDelays),trial_nr);
    seed = 8888;
    rng(seed)
    for i=1:trial_nr
%         if i == 1
%             pidx = seed;
%             qidx = seed+1;
%             f0idx = 1/seed;
%         end
%         rng(f0idx);
        f0 = round(2 * (rand(1, 1) - 0.5) * 5) + 65;
        go = true;
        while go
%             % broadband noise
%             rng(pidx); 
%             ya = randn(size(t));
%             rng(qidx); 
%             yb = randn(size(t));
%             
%             ya = filtfilt(bl, al, ya); 
%             ya = filtfilt(bh, ah, ya);
%             yb = filtfilt(bl, al, yb); 
%             yb = filtfilt(bh, ah, yb);
% 
%             tag1 = cos(angle(hilbert(ya)));
%             tag2 = cos(angle(hilbert(yb)));
            
            % phase-modulated signal
%             rng(pidx); 
            p = randn(size(t));
%             rng(qidx); 
            q = randn(size(t));
            
            p = filtfilt(bl, al, [p(end:-1:1), p, p(end:-1:1)]); p = p((L + 1):(2 * L)); %low
            p = filtfilt(bh, ah, [p(end:-1:1), p, p(end:-1:1)]); p = p((L + 1):(2 * L)); %high
            q = filtfilt(bl, al, [q(end:-1:1), q, q(end:-1:1)]); q = q((L + 1):(2 * L)); %low
            q = filtfilt(bh, ah, [q(end:-1:1), q, q(end:-1:1)]); q = q((L + 1):(2 * L)); %high
            
            p = p - min(p); p = p / max(p); p = 2 * pi * (p - 0.5);
            q = q - min(q); q = q / max(q); q = 2 * pi * (q - 0.5);
            
            tag1 = sin(2 * pi * t * f0 + p); %left hemifield
            tag2 = sin(2 * pi * t * f0 + q); %right hemifield
            
            % %% CONTROLS
            % tag1x = xcorr(tag1',100); 
            % figure(1); plot(abs(tag1x)); hold on; title(mean(abs(tag1x(abs(tag1x) > max(abs(tag1x))/2))));
            % tag2x = xcorr(tag2',100);
            % figure(2); plot(abs(tag2x)); hold on; title(mean(abs(tag2x(abs(tag2x) > max(abs(tag2x))/2))));
            % 
            % subplot(1, 2, 1); plot(t, tag1 - 1); hold on; plot(t, tag1 + 1); plot(t, 0.1 * p); plot(t, 0.1 * q); xlabel('t (s)'); xlim([0, 1.0]); title('frequency modulation');
            % title(sprintf('f0=%d', f0), 'FontWeight', 'normal', 'FontSize', 10);
            % subplot(1, 2, 2); plot(f, abs(fft(tag1))); hold on; plot(f, abs(fft(tag1))); xlim([0, 100]); xlabel('f (Hz)');
            % %%
            
%             pidx = pidx+2;
%             qidx = qidx+2;
%             f0idx = 1/pidx;
            
            tag1x = xcorr(tag1',100); 
            %plot(abs(tag1x)); title(sum(abs(tag1x(abs(tag1x) > max(abs(tag1x))/2)))/max(abs(tag1x)));
            if sum(abs(tag1x(abs(tag1x) > max(abs(tag1x))/2)))/max(abs(tag1x)) < 10
                tag2x = xcorr(tag2',100);
                %plot(abs(tag2x)); title(sum(abs(tag2x(abs(tag2x) > max(abs(tag2x))/2)))/max(abs(tag2x)));
                if sum(abs(tag2x(abs(tag2x) > max(abs(tag2x))/2)))/max(abs(tag2x)) < 10
                    if abs(corr(tag1(:), tag2(:))) < 0.01
                        tag_list{1,i} = tag1;
                        tag_list{2,i} = tag2;
                        %figure(1); plot(abs(tag1x)); hold on; title(mean(abs(tag1x(abs(tag1x) > max(abs(tag1x))/2))));
                        %figure(2); plot(abs(tag2x)); hold on; title(mean(abs(tag2x(abs(tag2x) > max(abs(tag2x))/2))));
                        go = false;
                    end
                end
            end
        end
    end
    
    % %% CONTROLS
    % for i=1:size(tag_list,2)
    %    tag1 = tag_list(1,i);
    %    tag2 = tag_list(2,i);
    %    tag_cor(i) = corr(tag1{:}', tag2{:}');    
    % end
    % 
    % tags1 = cell2mat(tag_list(1,:)');
    % tags2 = cell2mat(tag_list(2,:)');
    % f_alltrials = linspace(0, frate, frate*2);
    % for i=1:trial_nr
    %     t1 = abs(fft(tags1(i,:),frate*2));
    %     t2 = abs(fft(tags2(i,:),frate*2));
    %     figure(1);plot(f_alltrials,t1); xlim([0, 100]); xlabel('f (Hz)'); hold on
    %     figure(2);plot(f_alltrials,t2); xlim([0, 100]); xlabel('f (Hz)'); hold on
    % end
    % 
    % tagss1 = cell2mat(tag_list(1,:));
    % tagss2 = cell2mat(tag_list(2,:));
    % figure; plot(abs(xcorr(tagss1',100))); hold on; plot(abs(xcorr(tagss2',100))); plot(abs(xcorr(tagss1',tagss2',100)));
    % a = xcorr(tagss1', tagss1');
    % figure; plot(f_alltrials, abs(fft(tagss1,frate*2))); hold on; plot(f_alltrials, abs(fft(tagss2,frate*2))); xlim([0, 100]); xlabel('f (Hz)');
    % %%
    
    % normalize signal
    for i=1:numel(tag_list)
        tag_list{i} = tag_list{i} - min(tag_list{i});
        tag_list{i} = tag_list{i} / max(tag_list{i});
    end
    
    % set the number of dimensions based on the projector mode (b/w or color)
    if strcmp(exp_version, 'bw')
        dim = 3;
    elseif strcmp(exp_version, 'color')
        dim = 1;
    end
    
    % generate 4 tagging signals (one per quadrant)
    for i=1:numel(tag_list)
        tag_list{i} = reshape(tag_list{i}, 4, dim, []);
    end
    
end

% abortExp function
function abortExp(debug, el)%, eyelink)

% log(outputfile, ['END of Block ' int2str(cfg.block) ' -ABORTED-']);

%Return propixx to normal state
if ~debug
    Datapixx('SetPropixxDlpSequenceProgram', 0);
    Datapixx('RegWrRd');
    Datapixx('close');
end

%lower priority
if ~debug
    Priority(0);
end

%stop eyelink
if el.eyelink
    Eyelink('Message', 'END OF SESSION - ABORTED');
    eyelinkStop(el);
end

%close all files
fclose('all');

%close screen
sca

%throw warning due to prematurely aborted experiment
warning('Experiment aborted');

end
