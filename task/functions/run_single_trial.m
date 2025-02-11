function[dataStr, win] = run_single_trial(td, scr, const,visual,symbols_textures, total_score, block_score)

    % draw trial information on EyeLink operator screen
    Eyelink('command','draw_cross %d %d', scr.xCenter, scr.yCenter);

    % predefine time stamps
    tOn    = NaN;
    tSac   = NaN; 
    ex_fg = NaN;
    
    %
    soa_range = [0.2, 0.5];
    maxRT = 3;

    % random
    soa = rand(1)* (soa_range(2)-soa_range(1)) + soa_range(1);

    % draw fixation & placeholders
    Screen('DrawDots', scr.main, visual.fix_location , round(scr.ppd*0.2), scr.black,[], 4); % fixation
    %Screen('DrawDots', scr.main, visual.tar_locations , visual.tar_size+round(scr.ppd/2), scr.lightgrey ,[], 2); % placeholders
    Screen('FillOval', scr.main, scr.lightgrey, visual.disc_rect);
    % draw coins in 1 go
    if block_score > 0
        Screen('DrawTextures', scr.main,  visual.token_tex, [], visual.rect_coin(:,1:block_score));
    end
    if total_score > 0
        DrawFormattedText(scr.main, num2str(total_score), visual.score_location(1), visual.score_location(2), scr.black);
    end
    tFix = Screen('Flip', scr.main,0);
    Eyelink('message', 'EVENT_FixationDot');

    % fixation check
    while GetSecs < (tFix + soa - scr.fd)
        [x,y] = getCoord(scr, const); % get eye position data
        chkres = checkGazeAtPoint([x,y],[scr.xCenter, scr.yCenter], visual.fixCkRad);
        if ~chkres
            ex_fg = 1;
        end
    end
    
    % draw stimuli 
    Screen('DrawDots', scr.main, visual.fix_location , round(scr.ppd*0.2), scr.black,[], 4); % fixation
     %Screen('DrawDots', scr.main, visual.tar_locations , visual.tar_size+round(scr.ppd/2), scr.lightgrey ,[], 2); % placeholders
    Screen('FillOval', scr.main, scr.lightgrey, visual.disc_rect);
    % Screen('DrawTextures', scr.main,  symbols_textures(td.index), [], visual.tar_rect(td.index,:)');
    Screen('DrawTextures', scr.main,  symbols_textures(td.index), [], visual.tar_rect');
    % draw coins in 1 go
    if block_score > 0
        Screen('DrawTextures', scr.main,  visual.token_tex, [], visual.rect_coin(:,1:block_score));
    end
    if total_score > 0
        DrawFormattedText(scr.main, num2str(total_score), visual.score_location(1), visual.score_location(2), scr.black);
    end
    tOn = Screen('Flip', scr.main);
    Eyelink('message', 'EVENT_TargetOnset');

    % loop gaze control for response period
    while GetSecs < (tOn + maxRT)

        if isnan(tSac)
            [x,y] = getCoord(scr, const); % get eye position data
            chkres = checkGazeAtPoint([x,y],[scr.xCenter, scr.yCenter],visual.fixCkRad);
            % fprintf('%.2f\t%.2f\n',x,y);
            if chkres==0
                tSac = GetSecs;
                Eyelink('message', 'EVENT_Saccade1Started');
                break;
            end
        end
        
    end
    
    % check which target was chosen
    tar_loc_ordered = visual.tar_locations(:,td.index);
    tar_choice = NaN;

    if isnan(tSac)
        ex_fg = 2;
    else

        while  GetSecs < (tOn + maxRT)
            
            [x,y] = getCoord(scr, const); % get eye position data

            % if checkGazeAtPoint([x,y],  tar_loc_ordered(:,1), visual.tar_ecc - 3*scr.ppd)==1
            if checkGazeAtPoint([x,y],  tar_loc_ordered(:,1), visual.tar_ecc - 3*scr.ppd)==1
                tar_choice = 1;
                ex_fg = 0;
                break;

            elseif checkGazeAtPoint([x,y],  tar_loc_ordered(:,2), visual.tar_ecc - 3*scr.ppd)
                tar_choice = 2;
                ex_fg = 0;
                break;

            end
        end
    end
    
    % here redraw with different color based on selection?
      % draw stimuli 
    Screen('DrawDots', scr.main, visual.fix_location , round(scr.ppd*0.2), scr.black,[], 4); % fixation
    %Screen('DrawDots', scr.main, visual.tar_locations(:,tar_choice) , visual.tar_size+round(scr.ppd), scr.colchosen ,[], 2); % 
    Screen('FillOval', scr.main, scr.colchosen, visual.choice_rect(:, td.index(tar_choice))); 
    %Screen('DrawDots', scr.main, visual.tar_locations , visual.tar_size+round(scr.ppd/2), scr.lightgrey ,[], 2); % placeholders
    Screen('FillOval', scr.main, scr.lightgrey, visual.disc_rect);
    % Screen('DrawTextures', scr.main,  symbols_textures(td.index), [], visual.tar_rect(td.index,:)');
    Screen('DrawTextures', scr.main,  symbols_textures(td.index), [], visual.tar_rect');
    % draw coins in 1 go
    if block_score > 0
        Screen('DrawTextures', scr.main,  visual.token_tex, [], visual.rect_coin(:,1:block_score));
    end
   if total_score > 0
        DrawFormattedText(scr.main, num2str(total_score), visual.score_location(1), visual.score_location(2), scr.black);
    end
    tOn = Screen('Flip', scr.main);
    Eyelink('message', 'EVENT_ChoiceComplete');
    
    % sample reward
    if ex_fg == 0
        P = td.probs(td.index(tar_choice));
        win = 0;
        if rand(1) <= P
            win = 1;
        end
    else
        P = -1;
        win =-1;
    end
    
    % win animation
    
    % this uses bezier curves to animate the motion of the coins
    if win > 0
        
%         if win >1
%             path_length = 20;
%         else
%             path_length = 40;
%         end
        
        path_length = 30;

        for p_i = 1:win
            
            % draw path
            [x_path, y_path] = bezierCurve2(visual.tar_locations(1,td.index(tar_choice)), visual.tar_locations(2,td.index(tar_choice)), visual.x_coins(block_score + 1),visual.y_coins, path_length);
            
            for i = 1:length(x_path)
                Screen('DrawDots', scr.main, visual.fix_location , round(scr.ppd*0.2), scr.black,[], 4); % fixation
                Screen('FillOval', scr.main, scr.colchosen, visual.choice_rect(:,td.index(tar_choice))); 
                Screen('FillOval', scr.main, scr.lightgrey, visual.disc_rect);
                %Screen('DrawTextures', scr.main,  symbols_textures(td.index), [], visual.tar_rect(td.index,:)');
                Screen('DrawTextures', scr.main,  symbols_textures(td.index), [], visual.tar_rect');
                
                % draw coins in 1 go
                if block_score > 0
                    Screen('DrawTextures', scr.main,  visual.token_tex, [], visual.rect_coin(:,1:block_score));
                end
                if total_score > 0
                    DrawFormattedText(scr.main, num2str(total_score), visual.score_location(1), visual.score_location(2), scr.black);
                end
                
                Screen('DrawTexture', scr.main, visual.token_tex, [], CenterRectOnPoint([0,0, visual.coin_size, visual.coin_size], round(x_path(i)), round(y_path(i))));
                Screen('Flip', scr.main);
            end
            
            % update score
            block_score = block_score + win;
            
            Screen('DrawDots', scr.main, visual.fix_location , round(scr.ppd*0.2), scr.black,[], 4); % fixation
            Screen('FillOval', scr.main, scr.colchosen, visual.choice_rect(:,td.index(tar_choice)));
            Screen('FillOval', scr.main, scr.lightgrey, visual.disc_rect);
            Screen('DrawTextures', scr.main,  symbols_textures(td.index), [], visual.tar_rect');
            
            if block_score > 0
                Screen('DrawTextures', scr.main,  visual.token_tex, [], visual.rect_coin(:,1:block_score));
            end
            if total_score > 0
                DrawFormattedText(scr.main, num2str(total_score), visual.score_location(1), visual.score_location(2), scr.black);
            end
            
            Screen('Flip', scr.main);
            
        end
        
%         if save_images
%             imageArray = Screen('GetImage', scr.main);
%             imwrite(imageArray, sprintf('./task_screenshots/im%i.jpg',im_counter));
%             im_counter = im_counter +1;
%         end
        
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

            % collect trial information TO BE FIXED
            trialData = sprintf('%i\t%.2f\t%i\t%i\t%i\t%.2f\t%.2f',[tar_choice P win block_score total_score td.probs ]);

            % add symbols names
            trialData = sprintf('%s\t%s\t%s', trialData, td.symbols{1}, td.symbols{2});
            
            % timing
            timeData = sprintf('%.2f\t%.2f\t%.2f',[tOn tSac soa]);

            % other response data
            rt = sprintf('%.2f',tSac - tOn);

            % collect data for tab [14 x trialData, 6 x timeData, 1 x respData]
            data = sprintf('%s\t%s\t%s\t%s',trialData, timeData, rt);

    end

    
    WaitSecs(0.5);
    
    Screen('DrawDots', scr.main, visual.fix_location , round(scr.ppd*0.2), scr.black,[], 4); % fixation
    Screen('FillOval', scr.main, scr.lightgrey, visual.disc_rect);
    
    if block_score > 0
        Screen('DrawTextures', scr.main,  visual.token_tex, [], visual.rect_coin(:,1:block_score));
    end
    if total_score > 0
        DrawFormattedText(scr.main, num2str(total_score), visual.score_location(1), visual.score_location(2), scr.black);
    end
    
    Screen('Flip', scr.main);
    
    dataStr = sprintf('%s',data); % print data to string
    %if const.TEST; fprintf(1,sprintf('\n%s',dataStr));end
    
    
 