% anti-saccade task; analysis single subject
clear all

addpath('../functions');
addpath('./analysis_functions');


%% screen settings
scr.subDist = 80;   % subject distance (cm)
scr.width   = 570;  % monitor width (mm)

scr.xres = 1920;
scr.yres = 1080;
scr.xCenter = scr.xres/2;
scr.yCenter = scr.yres/2 ;
ppd = va2pix(1,scr); % pixel per degree

%% import file
% location of raw data file 2300 2301
raw_data = '../data/AA01/S12.edf';

% system('edf2asc ../data/S1.edf -s -miss -1.0')

% load eye movement file
ds = edfmex(raw_data); % ,'-miss -1.0'

% see the content of the data
ds.FSAMPLE
ds.FEVENT.message

% how many trials? here are the index of img onsets for each trial
find(strcmp({ds.FEVENT.message}, 'EVENT_TargetOnset')==1)


%% prepare data

% which eye was tracked?
% 0=left 1=right (add 1 for indexing below)
eye_tracked = 2; %1 + mode([ds.FEVENT.eye]);

% initialize values of trial variables to NaN
trial_n = NaN;
trial_n_2 = NaN;
trial_n_3 = NaN;
t_start = NaN;
t_end =  NaN;
fixation_onset = NaN;
target_onset =  NaN;
choice_complete_tms =  NaN;
Soa= NaN;
block_n = NaN;
id = NaN;
tar_choice = NaN;
P = NaN;
win = NaN;
block_score= NaN; 
total_score = NaN;
prob_1= NaN;
prob_2= NaN;
img_1 = NaN;
img_2= NaN;
timestamp  =  [];
eye_x =  [];
eye_y =  [];

ds2 = {};
trial_count = 0;
% loop over all events in the file and extract relevant informations
% i = length(ds.FEVENT)-3
for i = 1:length(ds.FEVENT)

    % check if the event contain a 'message'
    % that is something included in the experimental code
    % with the command:  Eyelink('message', '<message>');
    % (this is tipycally used to record the timing of events, e.g. stimulus onset)
    if ~isempty(ds.FEVENT(i).message)

        % the message is a string and we
        % parse it in different "words"
        sa = strread(ds.FEVENT(i).message,'%s');

        % onset of trial
        if strcmp(sa(1),'TRIAL_START')
            trial_n = str2double(sa(2));
            t_start = ds.FEVENT(i).sttime;
        end

        % image onset
        if strcmp(sa(1),'EVENT_FixationDot')
            fixation_onset = int32(ds.FEVENT(i).sttime);
        end

        % EVENT_TargetOnset
        if strcmp(sa(1),'EVENT_TargetOnset')
            target_onset = int32(ds.FEVENT(i).sttime);
        end
        
         % EVENT_ChoiceComplete
        if strcmp(sa(1),'EVENT_ChoiceComplete')
            choice_complete_tms = int32(ds.FEVENT(i).sttime);
        end

        % end of eye movement recoding
        if strcmp(sa(1),'TRIAL_END')
            trial_n_2 = str2double(sa(2));
            t_end = ds.FEVENT(i).sttime;
        end

        if length(sa)>3

            % data info
            if strcmp(sa(1),'TrialData')

                id = sa(2);
                block_n = str2double(sa(3));
                trial_n_3 = str2double(sa(4));
                Soa = target_onset - fixation_onset;
                tar_choice = str2double(sa(5));
                P = str2double(sa(6));
                win = str2double(sa(7));
                block_score= str2double(sa(8));
                total_score = str2double(sa(9));
                prob_1= str2double(sa(10));
                prob_2= str2double(sa(11));
                img_1 = sa(12);
                img_2= sa(13);

            end
        end

    end

    % if we have everything, then extract gaze position samples
    if ~isnan(trial_n) && ~isnan(trial_n_2) && ~isnan(trial_n_3) && ~isnan(target_onset) 

        trial_count = trial_count+1;

        % find onset-offset of relevant recording
        % I noticed that some ms were missing in a trial, so I check for a
        % time 10 ms before I sent the message trial end
        index_start = find(ds.FSAMPLE.time==t_start);
        index_end = find(ds.FSAMPLE.time==(t_end-10));

        % timestamp (set 0 for target onset)
        timestamp  =  int32(ds.FSAMPLE.time(index_start:index_end)) - target_onset;

        eye_x =  ds.FSAMPLE.gx(eye_tracked, index_start:index_end);
        eye_y =  ds.FSAMPLE.gy(eye_tracked, index_start:index_end);

        % remove missing
        eye_x(eye_x==100000000 | eye_y==100000000) = NaN;
        eye_y(eye_y==100000000 | eye_x==100000000) = NaN;

        % remove blinks etc
        eye_x(eye_x<0 | eye_x>scr.xres | eye_y<0 | eye_y>scr.yres ) = NaN;
        eye_y(eye_y<0 | eye_y>scr.yres | eye_x<0 | eye_x>scr.xres ) = NaN;

        % save everything into a single structure
        ds2.trial(trial_count).trial_n = trial_n;
        ds2.trial(trial_count).t_start = t_start;
        ds2.trial(trial_count).t_end = t_end;
        ds2.trial(trial_count).fixation_onset = fixation_onset;
        ds2.trial(trial_count).target_onset = target_onset;
        ds2.trial(trial_count).Soa =  Soa;
        ds2.trial(trial_count).id = id;
        ds2.trial(trial_count).block_n = block_n;
        ds2.trial(trial_count).Soa = Soa;
        ds2.trial(trial_count).tar_choice = tar_choice;
        ds2.trial(trial_count).P = P;
        ds2.trial(trial_count).win = win;
        ds2.trial(trial_count).block_score= block_score;
        ds2.trial(trial_count).total_score = total_score;
        ds2.trial(trial_count).prob_1 = prob_1;
        ds2.trial(trial_count).prob_2 = prob_2;
        ds2.trial(trial_count).img_1 = img_1;
        ds2.trial(trial_count).img_2 = img_2;
        ds2.trial(trial_count).timestamp  =  timestamp  ;
        ds2.trial(trial_count).eye_x =  eye_x;
        ds2.trial(trial_count).eye_y =  eye_y;

        % re-initialize
        trial_n = NaN;
        trial_n_2 = NaN;
        trial_n_3 = NaN;
        t_start = NaN;
        t_end =  NaN;
        fixation_onset = NaN;
        target_onset =  NaN;
        choice_complete_tms =  NaN;
        Soa= NaN;
        block_n = NaN;
        id = NaN;
        tar_choice = NaN;
        P = NaN;
        win = NaN;
        block_score= NaN;
        total_score = NaN;
        prob_1= NaN;
        prob_2= NaN;
        img_1 = NaN;
        img_2= NaN;
        timestamp  =  [];
        eye_x =  [];
        eye_y =  [];

    end

end





