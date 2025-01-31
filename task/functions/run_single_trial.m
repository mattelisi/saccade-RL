function[dataStr] = run_single_trial(td, scr, const, design)

    % draw trial information on EyeLink operator screen
    Eyelink('command','draw_cross %d %d', scr.xCenter, scr.yCenter);

    % predefine time stamps
    tOn    = NaN;
    tSac   = NaN; 

    % random
    soa = rand(1)* (soa_range(2)-soa_range(1)) + soa_range(1);

    % draw fixation & placeholders
    Screen('DrawDots', scr.main, fix_location , round(ppd*0.2), scr.black,[], 4); % fixation
    Screen('DrawDots', scr.main, tar_locations , tar_size+round(ppd/2), scr.lightgrey ,[], 2); % placeholders
    tFix = Screen('Flip', scr.main,0);
    Eyelink('message', 'EVENT_FixationDot');

    % fixation check
    while GetSecs < (tFix + soa - scr.fd)
        [x,y] = getCoord(scr, const); % get eye position data
        chkres = checkGazeAtPoint([x,y],[scr.xCenter, scr.yCenter],fixCkRad);
        if ~chkres
            ex_fg = 1;
        end
    end
    
    % draw stimuli 
    Screen('DrawDots', scr.main, fix_location , round(ppd*0.2), scr.black,[], 4); % fixation
    Screen('DrawDots', scr.main, tar_locations , tar_size+round(ppd/2), scr.lightgrey ,[], 2); % placeholders
    Screen('DrawTextures', scr.main,  symbols_textures(td.index), [], tar_rect(td.index,:)');
    tOn = Screen('Flip', scr.main);
    Eyelink('message', 'EVENT_TargetOnset');

    % loop gaze control for response period
    while GetSecs < (tOn + maxRT)

        if isnan(tSac)
            [x,y] = getCoord(scr, const); % get eye position data
            chkres = checkGazeAtPoint([x,y],[scr.xCenter, scr.yCenter],fixCkRad);
            % fprintf('%.2f\t%.2f\n',x,y);
            if chkres==0
                tSac = GetSecs;
                Eyelink('message', 'EVENT_Saccade1Started');
                break;
            end
        end
        
    end
    
    % check which target was chosen
    tar_loc_ordered = tar_locations(:,td.index);
    tar_choice = NaN;

    if isnan(tSac)
        ex_fg = 2;
    else

        while  GetSecs < (tSac + 0.05)
            
            [x,y] = getCoord(scr, const); % get eye position data

            if checkGazeAtPoint([x,y],  tar_loc_ordered(:,1), tar_ecc - 3*ppd)==1
                tar_choice = 1;
                break;

            elseif checkGazeAtPoint([x,y],  tar_loc_ordered(:,2), tar_ecc - 3*ppd)
                tar_choice = 2;
                break;

            end
        end
    end

    if correct ==0
        beep;
    end

    %% trial end

    switch ex_fg

        case 1
            data = 'fixBreak';
            Eyelink('command','draw_text 100 100 15 Fixation break');

        case 2
            data = 'tooSlow';
            Eyelink('command','draw_text 100 100 15 Too slow or no saccade');

        case 0

            % collect trial information
            trialData = sprintf('%.2f\t%.2f\t%.2f\t',[soa target_side anti_saccade correct]);

            % timing
            timeData = sprintf('%.2f\t%.2f\t%.2f',[tFix tOn tSac]);

            % other response data
            rt = sprintf('%.2f',tSac - tOn);

            % collect data for tab [14 x trialData, 6 x timeData, 1 x respData]
            data = sprintf('%s\t%s\t%s\t%s',trialData, timeData, rt);

    end

    Eyelink('message', 'TRIAL_END %d',  t);
    Eyelink('stoprecording');

    dataStr = sprintf('%i\t%s\n',t,data); % print data to string
    %if const.TEST; fprintf(1,sprintf('\n%s',dataStr));end