function [] = run_experiment()
%---------------------------------------------------------------------------------------
%
% saccade RL task
% Matteo Lisi, 2025
%
%----------------------------------------------------------------------------------------
% Screen setup info: change these accordingto the monitor and viewing distance used
scr.subDist = 80;   % subject distance (cm)
scr.width   = 570;  % monitor width (mm)

%----------------------------------------------------------------------------------------
%% design
design.n_trials = 20; % per block
design.n_blocks = 6;

%----------------------------------------------------------------------------------------
% focus on the command window
commandwindow;
home;

addpath('./functions');

% Setup PTB with some default values
% PsychDefaultSetup(2);

% Seed the random number generator.
rng('shuffle')

% Skip sync tests for demo purposes only
Screen('Preference', 'SkipSyncTests', 2);

%----------------------------------------------------------------------------------------
%% collect some info?
vpcode = getVpCode;
SJ = getSJinfo;
if SJ.number > 0
    % info_str = sprintf('S%i', SJ.number);
    info_str = sprintf('%sS%i', vpcode, SJ.number);
    filename = sprintf('%sS%i', vpcode, SJ.number);
    session_n = SJ.number;
end

% create data fid
datFid = fopen(filename, 'w');

%----------------------------------------------------------------------
%% Screen setup

% Set the screen number to the external secondary monitor if there is one
% connected
screenNumber = max(Screen('Screens'));

% get rid of PsychtoolBox Welcome screen
%Screen('Preference', 'VisualDebugLevel',3);

% Define black, white and grey
scr.white = WhiteIndex(screenNumber);
scr.black = BlackIndex(screenNumber);
scr.grey = round(scr.white/2);
scr.lightgrey = scr.grey + 25;
scr.colchosen = 190;

% Open the screen
%[scr.main, scr.rect] = PsychImaging('OpenWindow', screenNumber, scr.grey);
[scr.main, scr.rect] = PsychImaging('OpenWindow', screenNumber, scr.grey, [0 0 1920 1080], 32, 2);
%imagingMode = kPsychNeed32BPCFloat;
%[scr.main, scr.rect] = Screen('OpenWindow',screenNumber, [0.5 0.5 0.5],[],32,2,0,2,imagingMode);
%[scr.main, scr.rect] = Screen('OpenWindow',screenNumber, scr.grey);

% Flip to cleartext_size
Screen('FillRect',scr.main, scr.grey);
Screen('Flip', scr.main);

% Query the frame duration
scr.fd = Screen('GetFlipInterval', scr.main);

% Query the maximum priority level
MaxPriority(scr.main);

% Get the centre coordinate of the window
[scr.xCenter, scr.yCenter] = RectCenter(scr.rect);
[scr.xres, scr.yres] = Screen('WindowSize', scr.main); % heigth and width of screen [pix]

% Set the blend funciton for the screen
Screen('BlendFunction', scr.main, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

ppd = va2pix(1,scr);   % pixel per degree conversion factor

text_size = round(ppd);%round((ppd/56) * 21); % text size


% -------------------------------------------------------------------
%% init eyelink stuff

const.TEST = 0; % dummy eyelink mode
visual.bgColor = [scr.grey,scr.grey,scr.grey];	% background color when calibrating
visual.fgColor = scr.black; % foreground color when calibrating
visual.ppd = ppd;

% initialize eyelink-connection
[el, err]=initEyelink(filename,visual,const,scr);
if err==el.TERMINATE_KEY
    return
end

% determine recorded eye
if ~isfield(const,'recEye') && ~const.TEST
    %evt = Eyelink('newestfloatsample');
    %const.recEye = find(evt.gx ~= -32768);
    eye_used = Eyelink('EyeAvailable'); % get tracked eye 
    const.recEye = Eyelink('EyeAvailable');
end

% tracker calibration
if ~const.TEST
    calibresult = EyelinkDoTrackerSetup(el);
    if calibresult==el.TERMINATE_KEY
        return
    end
end

%----------------------------------------------------------------------
%% visual settings

tar_ecc = ppd*7;
fix_location = [scr.xCenter, scr.yCenter];
tarX_locations = round([scr.xCenter - tar_ecc, scr.xCenter + tar_ecc]);
tar_locations = [tarX_locations; scr.yCenter, scr.yCenter];
tar_size = round(4*ppd);

fixCkRad = round(2*ppd);

tar_rect(1,:)= CenterRectOnPoint([0,0, tar_size, tar_size], tar_locations(1,1), tar_locations(2,1));
tar_rect(2,:)= CenterRectOnPoint([0,0, tar_size, tar_size], tar_locations(1,2), tar_locations(2,2));


%% Time settings
tFix = 0.2;
soa_range = [0.2, 0.4];
maxRT = 5;

%% img/textures to add:
imageFolder = [pwd, '/img/']; % 

% Get list of all files in the folder
imageFiles = dir(fullfile(imageFolder, '*.png')); % Adjust if there are other formats

% Extract filenames
fileNames = {imageFiles.name};

% Find the token image
token_name = 'mariostarcoin1.png';
token_path = fullfile(imageFolder, token_name);

% Store all other image paths
% symbols_paths = fullfile(imageFolder, fileNames(~strcmp(fileNames, token_name)));
symbols_names = fileNames(~strcmp(fileNames, token_name));

% Shuffle symbols randomly
rng('shuffle'); % Ensure randomness across runs
shuffled_symbols = symbols_names(randperm(length(symbols_names)));

if mod(length(shuffled_symbols), 2) ~= 0
    shuffled_symbols = shuffled_symbols(1:length(shuffled_symbols)-1);
end

% Create pairs without replacement
pairs = reshape(shuffled_symbols, 2, [])'; % Each row is a unique pair

% Define probabilities
probabilities = [0.25, 0.75]; % One symbol will be 0.25, the other 0.75

% Initialize design structure
design.b = struct();

% Assign pairs to blocks
for b = 1:design.n_blocks

    % Get the pair for this block (cycling through pairs)
    pair_idx = mod(b-1, size(pairs, 1)) + 1;
    symbol_pair = pairs(pair_idx, :);

    % Randomly assign high (0.8) and low (0.2) probability to symbols
    if rand < 0.5
        block_probs = probabilities; % [0.2, 0.8]
    else
        block_probs = fliplr(probabilities); % [0.8, 0.2]
    end

    % Store block-level assignment
    design.b(b).symbols = symbol_pair;
    design.b(b).probs = block_probs;

    % Initialize trials for the block
    trials = struct('symbols', {}, 'probs', {}, 'index', {});

    % Generate trials with pseudorandom order
    for t = 1:design.n_trials
        if rand < 0.5
            trials(t).symbols = symbol_pair; % Original order
            trials(t).probs = block_probs;   % Match probabilities order
            trials(t).index = [1, 2];
        else
            trials(t).symbols = fliplr(symbol_pair); % Flipped order
            trials(t).probs = fliplr(block_probs);   % Flip probabilities too
            trials(t).index = [2, 1];
        end
    end

    % Shuffle trials within the block
    shuffled_indices = randperm(design.n_trials);
    design.b(b).t = trials(shuffled_indices);
end


%% coin locations
% these are for a single block, and then resetted in the next
% while the score adds up
coin_size = round(1*ppd);
max_coins = design.n_trials;
x_coins = round(linspace(1.5*coin_size, round(scr.xres/2-0.5*coin_size), max_coins));
y_coins = round(scr.yres - coin_size);
rect_coin = zeros(max_coins,4);
for i=1:max_coins
    rect_coin(i,:)= CenterRectOnPoint([0,0, coin_size, coin_size], x_coins(i), y_coins);
end

[token_img,~,alpha] = imread([pwd, '/img/' token_name]);
token_img(:, :, 4) = alpha;
token_tex = Screen('MakeTexture', scr.main, token_img); % make opengl texture out of image


%% make visual setting structure 
visual.tar_locations = tar_locations;
visual.tar_size = tar_size;
visual.fixCkRad = fixCkRad;
visual.tar_rect = tar_rect;
visual.fix_location = fix_location;
visual.tar_ecc = tar_ecc;

%
scr.ppd = ppd; 
visual.text_size = text_size;

visual.disc_rect = [CenterRectOnPoint([0,0, tar_size+round(scr.ppd/2), tar_size+round(scr.ppd/2)], tar_locations(1,1),  tar_locations(2,1))', ...
    CenterRectOnPoint([0,0, tar_size+round(scr.ppd/2), tar_size+round(scr.ppd/2)], tar_locations(1,2),  tar_locations(2,2))'];

visual.choice_rect = [CenterRectOnPoint([0,0, tar_size+round(scr.ppd), tar_size+round(scr.ppd)], tar_locations(1,1),  tar_locations(2,1))', ...
    CenterRectOnPoint([0,0, tar_size+round(scr.ppd), tar_size+round(scr.ppd)], tar_locations(1,2),  tar_locations(2,2))'];

% coins info
visual.x_coins = x_coins;
visual.y_coins = y_coins;
visual.rect_coin = rect_coin';
visual.token_tex = token_tex;
visual.score_location = round([scr.xres-2*coin_size, visual.y_coins + text_size/2]);
visual.coin_size = coin_size;


% ------------------------------------------------------------------
%% misc
% Make a directory for the results
resultsDir = [pwd '/data/'];
if exist(resultsDir, 'dir') < 1
    mkdir(resultsDir);
end

%----------------------------------------------------------------------
%% run experiment

% total score
total_score = 0;

% wait for a key press to start
instructions = 'Press any key to begin';
Screen('TextSize', scr.main, text_size);
Screen('TextFont', scr.main, 'Arial');
DrawFormattedText(scr.main, instructions, scr.xCenter - ceil(scr.xCenter/1.2), 'center', scr.white);
Screen('Flip', scr.main);
SitNWait;

for b = 1:design.n_blocks

    % load textures
    [symbol_1,~,alpha_1] = imread( char(fullfile(imageFolder,  design.b(b).symbols(1))));
    [symbol_2,~,alpha_2] = imread( char(fullfile(imageFolder,  design.b(b).symbols(2))));
    symbol_1(:, :, 4) = alpha_1;
    symbol_2(:, :, 4) = alpha_2;
    tex_1 = Screen('MakeTexture', scr.main, symbol_1); % make opengl texture
    tex_2 = Screen('MakeTexture', scr.main, symbol_2); 
    symbols_textures = [tex_1 ,tex_2];
    
    %
    block_score = 0;

    for t = 1:design.n_trials

        % This supplies a title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message '' Trial %d of %d''', t, design.n_trials - t);

        % this marks the start of the trial
        Eyelink('message', 'TRIALID %d', t);

        % get trial specification
        td = design.b(b).t(t);

        % cehck fixation before starting
        ncheck = 0;
        fix    = 0;
        record = 0;
        if const.TEST < 2
            while fix~=1 || ~record
                if ~record
                    Eyelink('startrecording');	% start recording
                    % You should always start recording 50-100 msec before required
                    % otherwise you may lose a few msec of data
                    WaitSecs(.1);
                    if ~const.TEST
                        key=1;
                        while key~= 0
                            key = EyelinkGetKey(el);		% dump any pending local keys
                        end
                    end

                    err=Eyelink('checkrecording'); 	% check recording status
                    if err==0
                        record = 1;
                        Eyelink('message', 'RECORD_START');
                    else
                        record = 0;	% results in repetition of fixation check
                        Eyelink('message', 'RECORD_FAILURE');
                    end
                end

                if fix~=1 && record

                    Eyelink('command','clear_screen 0');
                    Screen('FillRect',scr.main, scr.grey);
                    Screen('Flip', scr.main);
                    WaitSecs(0.1);

                    % CHECK FIXATION
                    fix = checkFix(scr, fixCkRad, const, fix_location, ppd);
                    ncheck = ncheck + 1;
                end

                if fix~=1 && record
                    % calibration, if maxCheck drift corrections did not succeed
                    if ~const.TEST
                        calibresult = EyelinkDoTrackerSetup(el);
                        if calibresult==el.TERMINATE_KEY
                            return
                        end
                    end
                    record = 0;
                end
            end
        else
            Screen('DrawDots', scr.main, fix_location , round(ppd*0.2), scr.white,[], 4); % fixation
            Screen('Flip', scr.main);
            WaitSecs(0.2);
        end

        Eyelink('message', 'TRIAL_START %d', t);
        Eyelink('message', 'SYNCTIME');		% zero-plot time for EDFVIEW


        %% run single trial here
        [dataStr, win]  = run_single_trial(td, scr, const,visual,symbols_textures, total_score, block_score);

        if win >0
            block_score= block_score+win;
        end
        
        dataline = sprintf('%s\t%i\t%i\t%i\t%s', info_str, session_n, b, t, dataStr);
        fprintf(datFid, dataline);
        
        % save trial info to eye mv rec
        Eyelink('message','TrialData %s', sprintf('%s\t%i\t%i\t%i\t%s', info_str, session_n, b, t, dataStr));        
        Eyelink('message', 'TRIAL_END %d',  t);
        Eyelink('stoprecording');

        % go to next trial if fixation was not broken
        if strcmp(dataStr,'fixBreak')
            trialDone = 0;
            feedback('Please maintain fixation until target appears.',tar_locations(1),tar_locations(2),scr,visual);
            
        elseif strcmp(dataStr,'tooSlow')
            trialDone = 0;
            feedback('Too slow.',tar_locations(1),tar_locations(2),scr,visual);
            
        else
            trialDone = 1;
            Eyelink('message', 'TrialData %s', dataStr);% write data to edfFile
            fprintf(datFid,dataStr);                    % write data to datFile
        end

        % isi
        WaitSecs(0.75);

    end
    
    
    % update score at the end of the block
    if b < design.n_blocks
    instructions = 'End of the block. The symbols will now change.\n\n Press a key to continue.';
    Screen('TextSize', scr.main, text_size);
    Screen('TextFont', scr.main, 'Arial');
    DrawFormattedText(scr.main, instructions, scr.xCenter - ceil(scr.xCenter/1.2), 'center', scr.white);
    if block_score > 0
        Screen('DrawTextures', scr.main,  visual.token_tex, [], visual.rect_coin(:,1:block_score));
    end
    if total_score > 0
        DrawFormattedText(scr.main, num2str(total_score), visual.score_location(1), visual.score_location(2), scr.black);
    end
    Screen('Flip', scr.main);
    
    WaitSecs(0.5);
    
    SitNWait;
    end
   
    b_count =  block_score; 
    path_length = 15;

    for b_i = 1:b_count
        
        % draw path
        [x_path, y_path] = bezierCurve2(visual.x_coins(block_score), visual.y_coins, visual.score_location(1),visual.score_location(2), path_length);
        
        for i = 1:length(x_path)
            
            DrawFormattedText(scr.main, instructions, scr.xCenter - ceil(scr.xCenter/1.2), 'center', scr.white);
            if block_score > 0
                Screen('DrawTextures', scr.main,  visual.token_tex, [], visual.rect_coin(:,1:block_score));
            end
            if total_score > 0
                DrawFormattedText(scr.main, num2str(total_score), visual.score_location(1), visual.score_location(2), scr.black);
            end
            
            Screen('DrawTexture', scr.main, token_tex, [], CenterRectOnPoint([0,0, coin_size, coin_size], round(x_path(i)), round(y_path(i))));
            Screen('Flip', scr.main);
        end
        
        % update score
        block_score = block_score - 1;
        total_score = total_score + 1;
        
        DrawFormattedText(scr.main, instructions, scr.xCenter - ceil(scr.xCenter/1.2), 'center', scr.white);
        if block_score > 0
            Screen('DrawTextures', scr.main,  visual.token_tex, [], visual.rect_coin(:,1:block_score));
        end
        if total_score > 0
            DrawFormattedText(scr.main, num2str(total_score), visual.score_location(1), visual.score_location(2), scr.black);
        end
        
        %Screen('DrawTexture', scr.main, token_tex, [], CenterRectOnPoint([0,0, coin_size, coin_size], round(x_path(i)), round(y_path(i))));
        Screen('Flip', scr.main);
        
        
    end
    
end

%----------------------------------------------------------------------
%% close data file
fclose(datFid); % close datFile

% save also mat file so we have everything
<<<<<<< HEAD
save(sprintf('%s.mat', info_str),'design','visual','scr','const');
=======
save(sprintf('S%i.mat',SJ.number),'design','visual','scr','const');
>>>>>>> a417f14a21de4603a760b78b71f3ce74cdb9b696

%----------------------------------------------------------------------
%% final feedback and end screen
Screen('TextSize', scr.main, text_size);
Screen('TextFont', scr.main, 'Arial');

DrawFormattedText(scr.main, 'The end! Thanks for participanting.\n\n\n(press any key to exit)', scr.xCenter - ceil(scr.xCenter/1.2), scr.yCenter, scr.white);
Screen('Flip', scr.main);
SitNWait;

% wrap up eyething stuff and receive data file
reddUp;

%----------------------------------------------------------------------
%% Close the onscreen window
Screen('CloseAll');
