%function batch_eyetracking_analysis()

    % -- Path settings (adjust as needed) --
    dataDir = '../data';          % Where data folders live
    participantListFile = 'PID_list';  % Text file with one participant ID per line
    outFileName = 'analysis_all_data.csv'; % Output CSV name
    
    % -- Read participant IDs from text file --
    fid = fopen(participantListFile,'r');
    pidList = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    pidList = pidList{1};  % cell array of participant IDs
    
    % -- Initialize a table to collect all data --
    allData = table();
    
    % -----------------------------
    % Screen settings and constants
    % -----------------------------
    scr.subDist  = 80;   % cm
    scr.width    = 570;  % mm
    scr.xres     = 1920;
    scr.yres     = 1080;
    scr.xCenter  = scr.xres/2;
    scr.yCenter  = scr.yres/2;
    % Convert 1° to pixels
    ppd = va2pix(1, scr);  % << ensure va2pix is in your path
    
    % Additional relevant variables
    tar_ecc     = ppd*7;
    fix_location = [scr.xCenter, scr.yCenter];
    tarX_locations = round([scr.xCenter - tar_ecc, scr.xCenter + tar_ecc]);
    tar_locations  = [tarX_locations; scr.yCenter, scr.yCenter];
    tar_size   = round(4*ppd);
    fixCkRad   = round(2*ppd);
    
    % Saccade detection parameters
    SAMPRATE  = 1000;    % Eyetracker sampling rate 
    velSD     = 5;       
    minDur    = 8;       % threshold duration for saccades (ms)
    VELTYPE   = 2;       
    maxMSAmp  = 1;       
    mergeInt  = 10;      
    
    % ------------------------------
    % Process each participant
    % ------------------------------
    for p = 1:numel(pidList)
        
        PID = pidList{p};  % e.g. 'AA01', 'HC01', etc.
        fprintf('Processing participant: %s\n', PID);
        
        % Locate the participant's folder
        participantDir = fullfile(dataDir, PID);
        if ~exist(participantDir, 'dir')
            warning('Folder does not exist: %s. Skipping.', participantDir);
            continue;
        end
        
        % Find all .edf files in this folder
        edfFiles = dir(fullfile(participantDir, '*.edf'));
        
        if isempty(edfFiles)
            fprintf('  No EDF files found for %s\n', PID);
            continue;
        end
        
        % For each .edf in this participant’s folder:
        for f = 1:numel(edfFiles)
            
            edfName = edfFiles(f).name;
            edfPath = fullfile(participantDir, edfName);
            fprintf('  Analyzing file: %s\n', edfName);
            
            % -- Load the EDF file, parse out relevant ds2 structure --
            ds = edfmex(edfPath);  % Adjust call as needed
            
            % (RE)Initialize structure to collect trial data from this file
            ds2 = struct();  
            ds2.trial = [];
            trial_count = 0;
            
            % Decide which eye was tracked
            % 0 = left, 1 = right, add +1 if needed.
            % Here we just set eye_tracked=2 as in your example. 
            eye_tracked = 2;  
            
            % Temporary storage for messages
            trial_n = NaN;  trial_n_2 = NaN;  trial_n_3 = NaN;
            t_start = NaN;  t_end    = NaN;
            fixation_onset = NaN;  target_onset = NaN; choice_complete_tms = NaN;
            Soa = NaN; block_n = NaN; id = NaN; tar_choice = NaN; 
            P = NaN; win = NaN; block_score = NaN; total_score = NaN; 
            prob_1 = NaN; prob_2 = NaN; img_1 = NaN; img_2 = NaN;
            timestamp = []; eye_x = []; eye_y = [];
            
            % -----------------------------
            % Extract trial-level data (ds2)
            % -----------------------------
            for i = 1:length(ds.FEVENT)
                
                if ~isempty(ds.FEVENT(i).message)
                    
                    sa = strsplit(ds.FEVENT(i).message);  % parse by space
                    % (Your original code used strread, you can adapt as needed)
                    
                    % 1) TRIAL_START
                    if strcmp(sa{1}, 'TRIAL_START')
                        trial_n = str2double(sa{2});
                        t_start = ds.FEVENT(i).sttime;
                    end
                    
                    % 2) EVENT_FixationDot
                    if strcmp(sa{1}, 'EVENT_FixationDot')
                        fixation_onset = int32(ds.FEVENT(i).sttime);
                    end
                    
                    % 3) EVENT_TargetOnset
                    if strcmp(sa{1}, 'EVENT_TargetOnset')
                        target_onset = int32(ds.FEVENT(i).sttime);
                    end
                    
                    % 4) EVENT_ChoiceComplete
                    if strcmp(sa{1}, 'EVENT_ChoiceComplete')
                        choice_complete_tms = int32(ds.FEVENT(i).sttime);
                    end
                    
                    % 5) TRIAL_END
                    if strcmp(sa{1}, 'TRIAL_END')
                        trial_n_2 = str2double(sa{2});
                        t_end = ds.FEVENT(i).sttime;
                    end
                    
                    % 6) TrialData ...
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
                    end
                end
                
                % If we have enough info to define a trial, extract the eye data
                if ~isnan(trial_n) && ~isnan(trial_n_2) && ~isnan(trial_n_3) && ~isnan(target_onset) && response==1
                    
                    trial_count = trial_count + 1;
                    
                    % find start-end indexes
                    idx_start = find(ds.FSAMPLE.time==t_start, 1);
                    idx_end   = find(ds.FSAMPLE.time==(t_end-10), 1);
                    
                    if isempty(idx_start) || isempty(idx_end)
                        % edge case: no valid indexes found
                        fprintf('    Warning: missing sample range for trial %d\n', trial_count);
                        
                        % re-init and skip
                        [trial_n, trial_n_2, trial_n_3, t_start, t_end, ...
                         fixation_onset, target_onset, choice_complete_tms, ...
                         Soa, block_n, id, tar_choice, P, win, ...
                         block_score, total_score, prob_1, prob_2, ...
                         img_1, img_2, timestamp, eye_x, eye_y] = deal(NaN);
                        continue;
                    end
                    
                    % Timestamps relative to target onset
                    timestamp = int32(ds.FSAMPLE.time(idx_start:idx_end)) - target_onset;
                    
                    % Gaze positions
                    gx = ds.FSAMPLE.gx(eye_tracked, idx_start:idx_end);
                    gy = ds.FSAMPLE.gy(eye_tracked, idx_start:idx_end);
                    
                    % Remove missing
                    gx(gx==100000000 | gy==100000000) = NaN;
                    gy(gy==100000000 | gx==100000000) = NaN;
                    
                    % Remove out-of-range
                    gx(gx < 0 | gx > scr.xres | gy < 0 | gy > scr.yres) = NaN;
                    gy(gy < 0 | gy > scr.yres | gx < 0 | gx > scr.xres) = NaN;
                    
                    % Store in ds2
                    ds2.trial(trial_count).trial_n       = trial_n;
                    ds2.trial(trial_count).t_start       = t_start;
                    ds2.trial(trial_count).t_end         = t_end;
                    ds2.trial(trial_count).fixation_onset= fixation_onset;
                    ds2.trial(trial_count).target_onset  = target_onset;
                    ds2.trial(trial_count).Soa           = Soa;
                    ds2.trial(trial_count).id            = id;
                    ds2.trial(trial_count).block_n       = block_n;
                    ds2.trial(trial_count).tar_choice    = tar_choice;
                    ds2.trial(trial_count).P             = P;
                    ds2.trial(trial_count).win           = win;
                    ds2.trial(trial_count).block_score   = block_score;
                    ds2.trial(trial_count).total_score   = total_score;
                    ds2.trial(trial_count).prob_1        = prob_1;
                    ds2.trial(trial_count).prob_2        = prob_2;
                    ds2.trial(trial_count).img_1         = img_1;
                    ds2.trial(trial_count).img_2         = img_2;
                    ds2.trial(trial_count).timestamp     = timestamp;
                    ds2.trial(trial_count).eye_x         = gx;
                    ds2.trial(trial_count).eye_y         = gy;
                    
                    % Re-initialize placeholders for the next trial
                    [trial_n, trial_n_2, trial_n_3, t_start, t_end, ...
                     fixation_onset, target_onset, choice_complete_tms, ...
                     Soa, block_n, id, tar_choice, P, win, ...
                     block_score, total_score, prob_1, prob_2, ...
                     img_1, img_2, timestamp, eye_x, eye_y] = deal(NaN);
                    
                end % end if we have a complete trial
            end % end loop over FEVENT
            
            
            % ---------------------------------------------------
            % Now ds2 has one entry per trial. Run the saccade analysis
            % ---------------------------------------------------
            
            for t = 1:numel(ds2.trial)
                
                % Extract the gaze data for this trial
                timers = ds2.trial(t).timestamp;
                XY = [ds2.trial(t).eye_x; ds2.trial(t).eye_y]';  % Nx2
                
                saccade_ok  = 0;
                sacOnset    = NaN;  sacOffset  = NaN;  sacDur     = NaN;
                sacVPeak    = NaN;  sacDist    = NaN;  sacAngle1  = NaN;
                sacAmp      = NaN;  sacAngle2  = NaN;  sacxOnset  = NaN;
                sacyOnset   = NaN;  sacxOffset = NaN;  sacyOffset = NaN;
                sacRT       = NaN;  sacChoice  = NaN;
                
                if ~isempty(timers)
                
                % Convert (X,Y) from pixels to degrees wrt screen center
                xrsf = double( (1/ppd) * [XY(:,1) - scr.xCenter, scr.yCenter - XY(:,2)] );
                
                % Filter (e.g., moving average)
                xrs = zeros(size(xrsf));
                xrs(:,1) = movmean(xrsf(:,1), 8);
                xrs(:,2) = movmean(xrsf(:,2), 8);
                
                % Compute velocities
                vrsf = vecvel(xrsf, SAMPRATE, VELTYPE);
                vrs  = vecvel(xrs,  SAMPRATE, VELTYPE);
                
                % Detect saccades
                mrs = microsaccMerge(xrsf, vrsf, velSD, minDur, mergeInt);
                mrs = saccpar(mrs);
                
                % Exclude micro-saccades smaller than maxMSAmp
                if ~isempty(mrs)
                    amp = mrs(:,7);
                    mrs = mrs(amp > maxMSAmp, :);
                end
                
                % Set up squares/circles for checking fixations
                fixRec = repmat([0,0],1,2)+[-2 -2 2 2];
                tarLRec = repmat([-7,0],1,2)+[-4 -4 4 4];
                tarRRec = repmat([ 7,0],1,2)+[-4 -4 4 4];
                
                
                
                % Find first saccade that leaves the fixation area and lands in L/R
                for s = 1:size(mrs,1)
                    xBeg  = xrs(mrs(s,1),1);
                    yBeg  = xrs(mrs(s,1),2);
                    xEnd  = xrs(mrs(s,2),1);
                    yEnd  = xrs(mrs(s,2),2);
                    
                    fixedFix = isincircle(xBeg,yBeg, fixRec);
                    landedL  = isincircle(xEnd,yEnd,  tarLRec);
                    landedR  = isincircle(xEnd,yEnd,  tarRRec);
                    
                    if fixedFix && (landedL || landedR)
                        saccade_ok = 1;
                        reaSacNumber = s;  %#ok<NASGU> 
                        
                        % Extract relevant saccade info
                        sacOnset   = mrs(s,1);
                        sacOffset  = mrs(s,2);
                        sacDur     = mrs(s,3);
                        sacVPeak   = mrs(s,4);
                        sacDist    = mrs(s,5);
                        sacAngle1  = mrs(s,6);
                        sacAmp     = mrs(s,7);
                        sacAngle2  = mrs(s,8);
                        
                        sacxOnset  = xrs(mrs(s,1),1);
                        sacyOnset  = xrs(mrs(s,1),2);
                        sacxOffset = xrs(mrs(s,2),1);
                        sacyOffset = xrs(mrs(s,2),2);
                        
                        % if needed, compute RT in ms: sacOnset index => timers(sacOnset)
                        %   because timers(1) = 0 means target onset
                        if (mrs(s,1) <= length(timers))
                            sacRT = timers(mrs(s,1));
                        end
                        
                        if landedR
                            sacChoice = 2;
                        elseif landedL
                            sacChoice = 1;
                        end
                        
                        break;  % first valid saccade
                    end
                end
                
                end
                
                % ----------------------------------------------
                % Collect a row for the final output table
                % ----------------------------------------------
                % Pull out DS2 trial-level fields we need:
                rowPID      = ds2.trial(t).id;            % from 'TrialData'  
                rowTrialN   = ds2.trial(t).trial_n;
                rowSoa      = ds2.trial(t).Soa;
                rowBlockN   = ds2.trial(t).block_n;
                rowTarChoice= ds2.trial(t).tar_choice;
                rowP        = ds2.trial(t).P;
                rowWin      = ds2.trial(t).win;
                rowBlockScr = ds2.trial(t).block_score;
                rowTotScr   = ds2.trial(t).total_score;
                rowProb1    = ds2.trial(t).prob_1;
                rowProb2    = ds2.trial(t).prob_2;
                rowImg1     = ds2.trial(t).img_1;
                rowImg2     = ds2.trial(t).img_2;
                
                % Build the table row
                T = table(...
                    {rowPID},         ... % keep as cell if it's a string
                    {edfName},        ...
                    rowTrialN,        ...
                    rowSoa,           ...
                    rowBlockN,        ...
                    rowTarChoice,     ...
                    rowP,             ...
                    rowWin,           ...
                    rowBlockScr,      ...
                    rowTotScr,        ...
                    rowProb1,         ...
                    rowProb2,         ...
                    {rowImg1},          ...
                    {rowImg2},          ...
                    sacOnset,         ...
                    sacOffset,        ...
                    sacDur,           ...
                    sacVPeak,         ...
                    sacDist,          ...
                    sacAngle1,        ...
                    sacAmp,           ...
                    sacAngle2,        ...
                    sacxOnset,        ...
                    sacyOnset,        ...
                    sacxOffset,       ...
                    sacyOffset,       ...
                    sacRT,            ...
                    sacChoice,        ...
                    'VariableNames', { ...
                        'participant_id','edf_file',   ...
                        'trial_n','Soa','block_n','tar_choice','P','win','block_score','total_score',...
                        'prob_1','prob_2','img_1','img_2',...
                        'sacOnset','sacOffset','sacDur','sacVPeak','sacDist','sacAngle1','sacAmp','sacAngle2',...
                        'sacxOnset','sacyOnset','sacxOffset','sacyOffset','sacRT','sacChoice'});
                
                % Append to master table
                allData = [allData; T];
                
            end % end trial loop
        end % end edf file loop
        
    end % end participant loop
    
    % -------------
    % Write CSV out
    % -------------
    writetable(allData, outFileName);
    fprintf('All done! Results saved to: %s\n', outFileName);

%end
