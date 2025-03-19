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

%% other task settings

tar_ecc = ppd*7;
fix_location = [scr.xCenter, scr.yCenter];
tarX_locations = round([scr.xCenter - tar_ecc, scr.xCenter + tar_ecc]);
tar_locations = [tarX_locations; scr.yCenter, scr.yCenter];
tar_size = round(4*ppd);

fixCkRad = round(2*ppd);

%% import file
% location of raw data file 2300 2301
raw_data = '../data/AA01/AA01S10.edf';

% system('edf2asc ../data/S1.edf -s -miss -1.0')

% load eye movement file
ds = edfmex(raw_data); % ,'-miss -1.0'

% see the content of the data
ds.FSAMPLE
ds.FEVENT.message

% % how many trials? here are the index of img onsets for each trial
% isx_allt = find(strcmp({ds.FEVENT.message}, 'TrialData')==1);
% isx_allt(61)

%% prepare data

% which eye was tracked?
% 0=left 1=right (add 1 for indexing below)
eye_tracked = 2; % 1 + mode([ds.FEVENT.eye]);

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

        if length(sa) > 3  && strcmp(sa{1}, 'TrialData') && strcmp(sa(5),'tooSlow')
            id        = sa{2};
            block_n   = str2double(sa{3});
            trial_n_3 = str2double(sa{4});
            Soa       = target_onset - fixation_onset;
            tar_choice = NaN;
            P         = NaN;
            win       = NaN;
            block_score = NaN;
            total_score = NaN;
            prob_1    = NaN;
            prob_2    = NaN;
            img_1     = NaN;
            img_2     = NaN;

            response = 0;

%             if (trial_count+1)==61 %debug
%                 print(sa);
%             end

        elseif length(sa) > 3 && strcmp(sa{1}, 'TrialData')  && ~strcmp(sa(5),'tooSlow')
            id        = sa{2};
            block_n   = str2double(sa{3});
            trial_n_3 = str2double(sa{4});
            Soa       = target_onset - fixation_onset;
            tar_choice = str2double(sa{5});
            P         = str2double(sa{6});
            win       = str2double(sa{7});
            block_score = str2double(sa{8});
            total_score = str2double(sa{9});
            prob_1    = str2double(sa{10});
            prob_2    = str2double(sa{11});
            img_1     = sa{12};
            img_2     = sa{13};

            response = 1;

%             if (trial_count+1)==61 %debug
%                 print(sa);
%             end

        end
    end

    % if we have everything, then extract gaze position samples
    if ~isnan(trial_n) && ~isnan(trial_n_2) && ~isnan(trial_n_3) && ~isnan(target_onset)  && response==1

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


%% saccade analysis for 1 image
t = 23;

% saccade algorithm parameters
SAMPRATE  = 1000;       % Eyetracker sampling rate 
velSD     = 5;          % lambda for microsaccade detectionc
minDur    = 8;          % threshold duration for microsaccades (ms)
VELTYPE   = 2;          % velocity type for saccade detection
maxMSAmp  = 1;          % maximum microsaccade amplitude
mergeInt  = 10;         % merge interval for subsequent saccadic events




xrs = [];
xrsf = [];
vrs = [];
vrsf= [];

% timestamps
timers = ds2.trial(t).timestamp;

% gaze position samples
XY = [ds2.trial(t).eye_x; ds2.trial(t).eye_y]';

% invert Y coordinates and transform in degrees relative to screen center
xrsf = double((1/ppd) * [XY(:,1)-scr.xCenter, (scr.yCenter-XY(:,2))]);  

% filter eye movement data (mostly for plotting)
xrs(:,1) = movmean(xrsf(:,1),8);
xrs(:,2) = movmean(xrsf(:,2),8);

% samrat = round(1000/mean(diff(ds2.trial(t).timestamp)));
% xrs(:,1) = filtfilt(fir1(35,0.05*SAMPRATE/samrat),1,xrsf(:,1));
% xrs(:,2) = filtfilt(fir1(35,0.05*SAMPRATE/samrat),1,xrsf(:,2));

% compute saccade parameters
vrs = vecvel(xrs, SAMPRATE, VELTYPE);    % velocities
vrsf= vecvel(xrsf, SAMPRATE, VELTYPE);   % velocities

mrs = microsaccMerge(xrsf,vrsf,velSD, minDur, mergeInt);  % saccades
mrs = saccpar(mrs);

% PLOT TRACES
% prepare figure and axes
close all;
cbac = [1.0 1.0 1.0];
h1 = figure;
set(gcf,'color',cbac);
ax(1) = axes('pos',[0.1 0.6 0.85 0.4]); % left bottom width height
ax(2) = axes('pos',[0.1 0.1 0.85 0.4]);

timers = double(ds2.trial(t).timestamp);
fix_onset = double(ds2.trial(t).fixation_onset - ds2.trial(t).target_onset)/1000;
% timeIndex = (timers - timers(1) +1)/1000;
timeIndex = (timers)/1000;


axes(ax(1));
% plot horizontal position
plot([0 0], [-12 12],'--','color',[0.5 0.5 0.5],'linewidth',0.8);
hold on
plot([fix_onset fix_onset], [-12 12],'--','color',[0.5 0.5 0.5],'linewidth',0.8);
plot(timeIndex,xrs(:,1),'-','color',[0.8 0 0],'linewidth',1);
for i = 1:size(mrs,1)
    plot(timeIndex((mrs(i,1):mrs(i,2))), xrs(mrs(i,1):mrs(i,2),1),'-','color',[0.8 0 0],'linewidth',3);
end
% plot vertical position
plot(timeIndex,xrs(:,2),'-','color',[0.2 0.2 0.8],'linewidth',1);
for i = 1:size(mrs,1)
    plot(timeIndex((mrs(i,1):mrs(i,2))), xrs(mrs(i,1):mrs(i,2),2),'-','color',[0 0 0.8],'linewidth',3);
end
ylim([-max(abs(xrs(:))),max(abs(xrs(:)))]*1.05)
ylabel('position [deg]');

axes(ax(2));
% plot horizontal vel
plot([0 0], [-600 600],'--','color',[0.5 0.5 0.5],'linewidth',1);
hold on
plot([fix_onset fix_onset], [-600 600],'--','color',[0.5 0.5 0.5],'linewidth',0.8);
plot(timeIndex,vrs(:,1),'-','color',[0.8 0 0],'linewidth',1);
for i = 1:size(mrs,1)
    plot(timeIndex((mrs(i,1):mrs(i,2))), vrs(mrs(i,1):mrs(i,2),1),'-','color',[0.8 0 0],'linewidth',3);
end
% plot vertical vel
plot(timeIndex,vrs(:,2),'-','color',[0.2 0.2 0.8],'linewidth',1);
for i = 1:size(mrs,1)
    plot(timeIndex((mrs(i,1):mrs(i,2))), vrs(mrs(i,1):mrs(i,2),2),'-','color',[0 0 0.8],'linewidth',3);
end
ylim([-max(abs(vrs(:))),max(abs(vrs(:)))]*1.05)

xlabel('time [sec]');
ylabel('velocity [deg/sec]');


% find the first saccade that leave fixation area 

% esclude microsasccades (amplitude smaller than cut off maxMSAmp)
if size(mrs,1)>0
    amp = mrs(:,7);
    mrs = mrs(amp>maxMSAmp,:);
end

% fixation check location
% fixRec_pix = repmat(fix_location,1,2)+[-fixCkRad -fixCkRad fixCkRad fixCkRad];
fixRec = repmat([0,0],1,2)+[-2 -2 2 2];
tarLRec = repmat([-7,0],1,2)+[-4 -4 4 4];
tarRRec = repmat([7,0],1,2)+[-4 -4 4 4];

saccade_ok = 0;

for s = 1:size(mrs,1)
    onset = timers(mrs(s,1));
    xBeg  = xrs(mrs(s,1),1);    % initial eye position x
    yBeg  = xrs(mrs(s,1),2);	% initial eye position y
    xEnd  = xrs(mrs(s,2),1);    % final eye position x
    yEnd  = xrs(mrs(s,2),2);	% final eye position y
    
    fixedFix = isincircle(xBeg,yBeg,fixRec);
    landedL = isincircle(xEnd,yEnd,tarLRec);
    landedR = isincircle(xEnd,yEnd,tarRRec);
    
    if fixedFix && (landedL||landedR) 
        saccade_ok = 1;  
        reaSacNumber = s;
        break;
    end
end

% other saccade parameters
sacOnset   = NaN;
sacOffset  = NaN;
sacDur     = NaN;
sacVPeak   = NaN;
sacDist    = NaN;
sacAngle1  = NaN;
sacAmp     = NaN;
sacAngle2  = NaN;
sacxOnset  = NaN;
sacyOnset  = NaN;
sacxOffset = NaN;
sacyOffset = NaN;
sacRT      = NaN;
sacChoice = NaN;

if saccade_ok == 1
    
    sacOnset   = mrs(reaSacNumber,1); % this is already relative to target onset
    sacOffset  = mrs(reaSacNumber,2);
    sacDur     = mrs(reaSacNumber,3); %*1000/samrat
    sacVPeak   = mrs(reaSacNumber,4);
    sacDist    = mrs(reaSacNumber,5);
    sacAngle1  = mrs(reaSacNumber,6);
    sacAmp     = mrs(reaSacNumber,7);
    sacAngle2  = mrs(reaSacNumber,8);
    sacxOnset  = xrs(mrs(reaSacNumber,1),1);
    sacyOnset  = xrs(mrs(reaSacNumber,1),2);
    sacxOffset = xrs(mrs(reaSacNumber,2),1);
    sacyOffset = xrs(mrs(reaSacNumber,2),2);
    
    if landedR == 1
        sacChoice = 2;
    elseif landedL == 1
        sacChoice = 1;
    end
end
