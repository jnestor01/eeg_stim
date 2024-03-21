function expt = rd_cupcakeAperture(subjectID, run)

%changes to make
%random seed (save in p)
%terminal output

%% Setup
if nargin==0
    subjectID = 'test';
    run = 1;
end

% % add paths
% addpath('../PTBWrapper/')
% addpath('../export-fig/')

HideCursor;

expName = 'CupcakeAperture';

saveData = 1;
plotTimingFigs = 1;
saveTimingFigs = 1;
plotStaircaseFigs = 1;
saveStaircaseFigs = 1;

p = cupcakeApertureParams;

% if strcmp(subjectID,'test')
%     p.eyeTracking = 0;
% end

if ~isempty(p.setSeed)
    rng(p.setSeed);
end

slack = (1/p.refRate)*.9;
ppd = ang2pix(1, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'); % pixels per degree

% Running on PTB-3? Abort otherwise.
AssertOpenGL;

if strcmp(p.testingLocation, 'desk')
    Screen('Preference', 'SkipSyncTests', 1);
end

%% Initialize stim tracker for MEG
if p.triggersOn
    PTBInitStimTracker;
    global PTBTriggerLength
    PTBTriggerLength = 0.001;
    
    % to send trigger messages to EyeLink
    global PTBEyeTrackerRecording
    if p.eyeTracking
        PTBEyeTrackerRecording = 1;
    else
        PTBEyeTrackerRecording = 0;
    end
end

%% Initialize parallel port triggers for EEG
if p.triggerEEG
    ppdev_mex('CloseAll');
    ppdev_mex('Open',1);
end

%% Display key settings to the experimenter
fprintf('\nExperiment settings:\n')
fprintf('aperture = %s\n', p.aperture)

% ok = input('\nSettings ok? [n if not]','s');
% if strcmp(ok,'n')
%     error('Ok, check parameters')
% end

%% Eye data i/o
eyeDataDir = 'eyedata';
eyeFile = sprintf('%s%02d%s', subjectID(1:2), run, datestr(now, 'mmdd'));
eyeFileFull = sprintf('%s/%s_run%02d_%s_%s.edf', eyeDataDir, subjectID, run, expName, datestr(now, 'yyyymmdd'));

% Check to see if this eye file already exists
existingEyeFile = dir(sprintf('%s/%s.edf', eyeDataDir, eyeFile));
if ~isempty(existingEyeFile) && p.eyeTracking
    error('eye file already exists! please choose another name.')
end

%% Check for existing data file
dataDir = 'data';
dataFile = sprintf('%s/%s_run%02d_%s_%s.mat', dataDir, subjectID, run, expName, datestr(now, 'yyyymmdd'));
existingDataFile = dir(dataFile);
if ~isempty(existingDataFile) && ~strcmp(subjectID, 'test') && ~strcmp(subjectID, 'testy')
    error('data file already exists!')
end

%% Keyboard
% set button box device
if p.useKbQueue
    devNum = [];
    devices = PsychHID('devices');
    for iD=1:numel(devices)
        if strcmp(devices(iD).usageName,'Keyboard') && ...
                strcmp(devices(iD).product,'904')
            devNum = iD;
        end
    end
    if isempty(devNum)
        error('Did not find button box')
    end
else
    devNum = -1;
end

% set up KbQueue if desired
if p.useKbQueue
    KbQueueCreate(devNum);
    KbQueueStart();
end

%% Sound
% Perform basic initialization of the sound driver
InitializePsychSound(1); % 1 for precise timing
PsychPortAudio('Close');
% Open audio device for low-latency output
pahandle = PsychPortAudio('Open', [], 1, 1, p.Fs, 2); % 1 = single-channel
Snd('Open', pahandle, 1); %nec for eyetracker


%% Screen
% Set up screen
screenNumber = max(Screen('Screens'));

% Check screen resolution and refresh rate
scr = Screen('Resolution', screenNumber);
% if ~all([scr.width scr.height scr.hz] == [p.screenRes p.refRate]) && ~strcmp(subjectID,'test')
%     error('Screen resolution and/or refresh rate has not been set correctly by the experimenter!')
% end

% Set up window
multisample = 8;
[window, rect] = Screen('OpenWindow', screenNumber,[],[],[],[],[],multisample);

% [window, rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]);
white = WhiteIndex(window);  % Retrieves the CLUT color code for white.
[cx, cy] = RectCenter(rect);
Screen('TextSize', window, p.fontSize);
Screen('TextColor', window, white);
Screen('TextFont', window, p.font);

% Check screen size
[sw, sh] = Screen('WindowSize', window); % height and width of screen (px)
% if ~all([sw sh] == p.screenRes) && ~strcmp(subjectID,'test')
%     error('Screen resolution is different from requested!')
% end

% Check refresh rate
flipInterval = Screen('GetFlipInterval', window); % frame duration (s)
% if abs(flipInterval - (1/p.refRate)) > 0.001 && ~strcmp(subjectID,'test')
%     error('Refresh rate is different from requested!')
% end

% Check font
if ~strcmp(p.font, Screen('TextFont', window))
    error('Font was not set to requested: %s', p.font)
end

% Load calibration file
switch p.testingLocation
    case 'CarrascoL1'
        load('../../Displays/Carrasco_L1_SonyGDM5402_sRGB_calibration_02292016.mat');
        table = CLUT;
        Screen('LoadNormalizedGammaTable', window, table);
        % check gamma table
        gammatable = Screen('ReadNormalizedGammaTable', window);
        if nnz(abs(gammatable-table)>0.0001)
            error('Gamma table not loaded correctly! Perhaps set screen res and retry.')
        end
    case 'MEG'
        d = load(sprintf('%s/gamma.mat', p.displayPath));
        table = d.gamma;
        Screen('LoadNormalizedGammaTable', window, table);
        % check gamma table
        gammatable = Screen('ReadNormalizedGammaTable', window);
        if nnz(abs(gammatable-table)>0.0001)
            error('Gamma table not loaded correctly! Perhaps set screen res and retry.')
        end
    otherwise
        fprintf('\nNot loading gamma table ...\n')
end

%% Load previous staircase value
prevDataFile = sprintf('%s/%s_run*_%s_*.mat', dataDir, subjectID, expName);
prevDataFiles = dir(prevDataFile);

if ~isempty(prevDataFiles)
    % find most recent file
    datenums = [];
    for i = 1:numel(prevDataFiles)
        datenums(i) = prevDataFiles(i).datenum;
    end
    [y, idx] = max(datenums);
    
    stairFileName = prevDataFiles(idx).name;
    stairFile = load(sprintf('%s/%s', dataDir, stairFileName));
    
    stairValues = stairFile.expt.staircase.stairValues;
    if isempty(stairValues) % there is a saved file, but the staircase was off
        stairIdx = numel(p.stairs); % start easy
        fprintf('\nStaircase start at %1.2f\n\n', p.stairs(stairIdx));
    else
        stairIdx = stairValues(end); % pick up where we left off
        fprintf('\nStaircase start at %1.2f, from file %s\n\n', p.stairs(stairIdx), stairFileName);
    end
    
    clear stairFile
else
    stairIdx = numel(p.stairs); % start easy
    fprintf('\nStaircase start at %1.2f\n\n', p.stairs(stairIdx));
end


%% Make stimuli
% Calculate stimulus dimensions (px) and position
imPos = round(p.imPos*ppd);
gratingRadius = round(p.gratingDiameter/2*ppd);
edgeWidth = round(p.apertureEdgeWidth*ppd);

% Make grating images and textures
orientation = 0;
spatialFrequency = p.gratingSF;

% p.imSize = repmat(p.gratingDiameter*2+1,[1,2]);
for iP = 1:numel(p.gratingPhases)
    phase = p.gratingPhases(iP);
    for iC = 1:numel(p.gratingContrasts)
        contrast = p.gratingContrasts(iC);
        
        grating = rd_grating(ppd, p.imSize(1), spatialFrequency, orientation, phase, contrast);
        
        ims{iP,iC} = rd_aperture(grating, p.aperture, gratingRadius, edgeWidth, p.angularFreq);
        texs(iP,iC) = Screen('MakeTexture', window, cat(3,ims{iP,iC}*white, (ims{iP,iC}~=0.5)*255));
    end
end

% Make the rects for placing the images
imSize = size(grating);

locs = round([cos(p.spatialLocations/180*pi)*ppd*p.radialDistance; sin(p.spatialLocations/180*pi)*ppd*p.radialDistance]);
imRect = repmat(CenterRectOnPoint([0 0 imSize], cx+imPos(1), cy+imPos(2)), numel(p.spatialLocations), 1);
imRect = imRect + [locs' locs'];

if p.showPlaceholders==1
    phRect = [min(imRect(:, [1 2]), [], 1), max(imRect(:,[3 4]), [], 1)] + [-1 -1 1 1]*p.phLineWidth;   
    phRect = phRect + [-1 -1 1 1]*round(p.phCushion*ppd); % expand placeholders by this many pixels so that they are not obscured by image rotations
elseif p.showPlaceholders==2
    scaleInner = (ppd*p.radialDistance-gratingRadius(1)-ppd*p.phCushion)/(ppd*p.radialDistance);
    scaleOuter = (ppd*p.radialDistance+gratingRadius(1)+ppd*p.phCushion)/(ppd*p.radialDistance);
    phRectInner = CenterRectOnPoint([round(min(locs, [], 2)'*scaleInner), max(locs, [], 2)'*scaleInner], cx, cy) + [-1 -1 1 1]*p.phLineWidth; 
    phRectOuter = CenterRectOnPoint([round(min(locs, [], 2)'*scaleOuter), max(locs, [], 2)'*scaleOuter], cx, cy) + [-1 -1 1 1]*p.phLineWidth;
    phRect = [phRectOuter; phRectInner];
end

% Where to put instruction text relative to center. For circular
% placeholders, it will be above the outer one or inside the inner one
% depending on where it fits
if p.showPlaceholders==2
    textRect = Screen('TextBounds', window, 'Experimenter press 5 to begin');
    textLoc = CenterRectOnPoint(textRect,cx,cy)+[0,round(3*ppd),0,round(3*ppd)];
    testRect(1:2) = textLoc(1:2)-[cx,cy];
    if sqrt(testRect(1)^2+testRect(2)^2)>ppd*(p.radialDistance*scaleInner-0.5)
        textY = cy-round(ppd*(p.radialDistance*scaleOuter+0.75));
    else
        textY = cy-round(3*ppd);
    end
else
    textY = cy-round(3*ppd);
end

% Calculate fixation diameter in pixels
p.fixSize = round(p.fixDiameter*ppd);

%% Generate trials
% Construct trials matrix
trials_headers = {'gratingOrientation','gratingPhase','gratingContrast','spatialLocation'...
    'targetState','iti','iti2','fixColor','responseKey','response','correct','rt'};

% make sure column indices match trials headers
gratingOrientationIdx = strcmp(trials_headers,'gratingOrientation');
gratingPhaseIdx = strcmp(trials_headers,'gratingPhase');
gratingContrastIdx = strcmp(trials_headers,'gratingContrast');
targetStateIdx = strcmp(trials_headers,'targetState');
itiIdx = strcmp(trials_headers,'iti');
iti2Idx = strcmp(trials_headers,'iti2');
fixColorIdx = strcmp(trials_headers,'fixColor');
responseKeyIdx = strcmp(trials_headers,'responseKey');
responseIdx = strcmp(trials_headers,'response');
correctIdx = strcmp(trials_headers,'correct');
rtIdx = strcmp(trials_headers,'rt');
spatialLocationIdx = strcmp(trials_headers,'spatialLocation');

% full factorial design
trials = fullfact([numel(p.gratingOrientations) ...
    numel(p.gratingPhases) ...
    numel(p.gratingContrasts) ...
    numel(p.spatialLocations)]);

% repeat trials matrix according to nReps of all conditions
trials = repmat(trials, p.nReps, 1);
nTrials = size(trials,1);

% show trials and blocks
fprintf('\n%s\n\n%d trials, %1.2f blocks\n\n', datestr(now), nTrials, nTrials/p.nTrialsPerBlock)

% Generate ITIs
switch p.itiType
    case 'uniform'
        nITIs = numel(p.itis);
        itis = repmat(1:nITIs,1,ceil(nTrials/nITIs));
        trials(:,itiIdx) = itis(randperm(nTrials));
    case 'hazard' % constant hazard rate
        trials(:,itiIdx) = rd_sampleDiscretePDF(p.itiPDF, nTrials);
    otherwise
        error('p.itiType not recognized')
end

% Generate target states
tsConds = ones(1,nTrials)*2; % target absent
idx = randperm(nTrials);
nTargets = round(nTrials*p.propTargetPresent);
tsConds(idx(1:nTargets)) = 1; % target present
trials(:,targetStateIdx) = tsConds;

%% Choose order of trial presentation
% trialOrder = randperm(nTrials);
% This randomizes trial order within reps, but not across reps. So we will
% present every trial in one rep before moving on to the next rep. In this
% way the experiment is blocked into complete reps (though the number of
% trials per block may be different then the number of trials in a rep).
nt = nTrials/p.nReps;
trialSet = ones(nt,1)*(1:p.nReps);
repOrders = [];
for i=1:p.nReps
    repOrders = [repOrders randperm(nt)'];
end
trialOrder0 = repOrders + (trialSet-1).*nt;
trialOrder = trialOrder0(:);

% shuffle trials
trials = trials(trialOrder,:);

%% Eyetracker
rad = round(ang2pix(p.eyeRad, p.screenSize(1), p.screenRes(1), p.viewDist, 'central'));
fixRect = [cx-rad, cy-rad, cx+rad, cy+rad];
if p.eyeTracking
    % Initialize eye tracker
    [el exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
    if exitFlag
        return
    end
    
    % Write subject ID into the edf file
    if ~strcmp(subjectID,'test')
    Eyelink('message', 'BEGIN DESCRIPTIONS');
    Eyelink('message', 'Subject code: %s', subjectID);
    Eyelink('message', 'Run: %d', run);
    Eyelink('message', 'END DESCRIPTIONS');
    end
    
    % No sounds indicating success of calibration
    el.drift_correction_target_beep = [0 0 0];
    el.drift_correction_failed_beep = [0 0 0];
    el.drift_correction_success_beep = [0 0 0];
    
    % Accept input from all keyboards
    el.devicenumber = -1; %see KbCheck for details of this value
    
    % Update with custom settings
    EyelinkUpdateDefaults(el);
    
    % Calibrate eye tracker
    [cal exitFlag] = rd_eyeLink('calibrate', window, el);
    if exitFlag
        return
    end
    rd_eyeLink('trialstart', window, {el, 1});

    % Start recording
    rd_eyeLink('startrecording',window,{el,fixRect});
end


Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
%% Example image and tones
% Wait for a button press

Screen('FillRect', window, p.backgroundColor*white);
drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders, p.ringColor*white);
drawFixation(window, cx, cy, p.fixSize, p.fixColor*white);

DrawFormattedText(window, 'Press 1 to start', 'center', textY, [1 1 1]*white);
Screen('Flip', window);

keyPressed = 0;
while ~keyPressed
    if p.useKbQueue
        [keyIsDown, firstPress] = KbQueueCheck();
        keyCode = logical(firstPress);
    else
        [secs, keyCode] = KbWait(devNum);
    end
    
    if strcmp(KbName(keyCode),'1!')
        keyPressed = 1;
    end
end

% Play example feedback tones
for iTone = 1:size(p.tones,1)
    PsychPortAudio('FillBuffer', pahandle, [p.tones(iTone,:); p.tones(iTone,:);]*p.soundAmp);
    PsychPortAudio('Start', pahandle, [], [], 1);
    WaitSecs(.8);
end

% Present example image
destRect = imRect(1,:);
tex = texs(1,2);
orientation = 0;
drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders, p.ringColor*white)
Screen('DrawTexture', window, tex, [], destRect, orientation);
drawFixation(window, cx, cy, p.fixSize, p.fixColor*white);

timeSampleIm = Screen('Flip', window);

% Show fixation
drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders, p.ringColor*white)
drawFixation(window, cx, cy, p.fixSize, p.fixColor*white);
Screen('Flip', window, timeSampleIm + p.gratingDur - slack);
WaitSecs(1);

% Now wait for the experimenter to initiate the run
% Show fixation
drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders, p.ringColor*white)
drawFixation(window, cx, cy, p.fixSize, p.fixColor*white);
DrawFormattedText(window, 'Experimenter press 5 to begin', 'center', textY, [1 1 1]*white);
Screen('Flip', window);

keyPressed = 0;
while ~keyPressed
    [secs, keyCode] = KbWait(-1);
    if strcmp(KbName(keyCode),'5%')
        keyPressed = 1;
    end
end

% Show fixation
drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders, p.ringColor*white)
drawFixation(window, cx, cy, p.fixSize, p.fixColor*white);
Screen('Flip', window, timeSampleIm + p.gratingDur - slack);
WaitSecs(2);

%% Present trials
% Trials
trialsPresented = [];
lastFewAcc = [];
stairValues = [];
reversalValues = [];
threshold = [];
triggerTimes = [];
extraITI = 0;
block = 1;
stairCounter = 1;

% Present initial fixation
drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders, p.ringColor*white)
drawFixation(window, cx, cy, p.fixSize, p.fixColor*white);
if strcmp(p.testingLocation,'MEG') || strcmp(p.testingLocation,'denisonlab-EEG')
    drawPhotodiode(window, [cx cy]*2, white, 0); % black
end
if p.triggerEEG
    ppdev_mex('Write', 1, p.triggers.fixation); 
end
timeFix = Screen('Flip', window); % first image will be presented after an ITI with respect to this time
timing.startTime = timeFix;
if p.triggerEEG
    ppdev_mex('Write', 1, 0); 
    triggerTimes = [triggerTimes; p.triggers.fixation GetSecs]; 
end
if p.triggersOn
    PTBSendTrigger(p.triggers.fixation, 0); 
    triggerTimes = [triggerTimes; p.triggers.fixation GetSecs];
elseif p.triggerEEG
%     ppdev_mex('Write', 1, p.triggers.fixation);
%     WaitSecs(0.005);
%     ppdev_mex('Write', 1, 0);
%     triggerTimes = [triggerTimes; p.triggers.fixation GetSecs];
end
if p.eyeTracking
    Eyelink('Message', 'EVENT_FIX');
end

% Present trials
for iTrial = 1:nTrials
    extraKeyPresses = [];
    if p.eyeTracking
        Eyelink('Message', 'TRIAL_START');
    end
    
    % Get conditions for this trial
    oriCond = trials(iTrial, gratingOrientationIdx);
    phaseCond = trials(iTrial, gratingPhaseIdx);
    contrastCond = trials(iTrial, gratingContrastIdx);
    tsCond = trials(iTrial, targetStateIdx);
    itiCond = trials(iTrial, itiIdx);
    spatialLocationCond = trials(iTrial, spatialLocationIdx);
    
    % Get actual values
    orientation = p.gratingOrientations(oriCond);
    phase = p.gratingPhases(phaseCond);
    contrast = p.gratingContrasts(contrastCond);
    targetState = p.targetStates(tsCond);
    iti = p.itis(itiCond) + extraITI;
    spatialLocation = p.spatialLocations(spatialLocationCond);
    
    % Select target textures
    tex = texs(phaseCond,contrastCond);

    % Select destination for given spatial location
    destRect = imRect(spatialLocationCond, :);
    
    % Set fixation contrast based on staircase
    if targetState==1 % target present
        if p.staircase
            fixColor = p.fixColor;
            fixColor(1) = p.stairs(stairIdx);
        else
            fixColor = p.stairs(end);
        end
    else
        fixColor = p.fixColor;
    end
    
    isFix = 1;
    eyeSlack = 0;
    frameTime = GetSecs;
    % Check keyboard during ITI
    while GetSecs < timeFix + iti - 3*slack + eyeSlack || ~isFix
        if p.useKbQueue
            [keyIsDown, firstPress] = KbQueueCheck();
            if keyIsDown
                secs = min(firstPress(firstPress~=0));
                keyCode = firstPress==secs;
            else
                keyCode = [];
            end
        else
            [keyIsDown, secs, keyCode] = KbCheck(devNum);
        end
        if keyIsDown
            rtx = secs - timeFix; % measure from start of trial
            key = find(keyCode);
            extraKeyPresses = [extraKeyPresses; secs rtx key(1)];
            
            if p.triggersOn
                PTBSendTrigger(p.triggers.response, 0);
                triggerTimes = [triggerTimes; p.triggers.response GetSecs];
            elseif p.triggerEEG
                ppdev_mex('Write', 1, p.triggers.response);
                WaitSecs(0.005);
                ppdev_mex('Write', 1, 0);
                triggerTimes = [triggerTimes; p.triggers.response GetSecs];
            end
        end
        if p.monitorFix
            isFix = rd_eyeLink('fixcheck', window, {cx, cy, rad});
            if ~isFix
                iti = iti+(GetSecs-frameTime);
                remain = timeFix + iti - 3*slack + eyeSlack - GetSecs;
                if remain < p.eyeSlack
                    eyeSlack = p.eyeSlack-remain;
                end
            end
        end
        frameTime = GetSecs;
    end
    iti = iti + eyeSlack;
    
%     Check fixation hold
%     if p.eyeTracking
%         fixation = 0;
%         while ~fixation
%             WaitSecs(.01);
%             fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
%         end
%     end
%     WaitSecs(p.eyeSlack);

    % Present image
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders, p.ringColor*white)
    Screen('DrawTexture', window, tex, [], destRect, orientation);
    drawFixation(window, cx, cy, p.fixSize, fixColor*white); % fixColor can change trial to trial
    if strcmp(p.testingLocation,'MEG') || strcmp(p.testingLocation,'denisonlab-EEG')
        drawPhotodiode(window, [cx cy]*2, white, 1); % white
    end
    if p.triggerEEG
        ppdev_mex('Write', 1, p.triggers.image); 
    end
    timeIm = Screen('Flip', window, timeFix + iti - slack);
    if p.triggerEEG 
        ppdev_mex('Write', 1, 0); 
        triggerTimes = [triggerTimes; p.triggers.image GetSecs]; 
    end
    if p.triggersOn
        PTBSendTrigger(p.triggers.image, 0);
        triggerTimes = [triggerTimes; p.triggers.image GetSecs];
        
        if targetState==1
            PTBSendTrigger(p.triggers.target, 0);
            triggerTimes = [triggerTimes; p.triggers.target GetSecs];
        end
    elseif p.triggerEEG
%         ppdev_mex('Write', 1, p.triggers.image);
%         WaitSecs(0.005);
%         ppdev_mex('Write', 1, 0);
%         triggerTimes = [triggerTimes; p.triggers.image GetSecs];

        if targetState==1
            ppdev_mex('Write', 1, p.triggers.target);
            WaitSecs(0.001);
            ppdev_mex('Write', 1, 0);
            triggerTimes = [triggerTimes; p.triggers.target GetSecs];
        end
    end
    if p.eyeTracking
        Eyelink('Message', 'EVENT_IMAGE');
    end
    
    % Check keyboard during image
    while GetSecs < timeIm + p.gratingDur - 3*slack
        if p.useKbQueue
            [keyIsDown, firstPress] = KbQueueCheck();
            if keyIsDown
                secs = min(firstPress(firstPress~=0));
                keyCode = firstPress==secs;
            else
                keyCode = [];
            end
        else
            [keyIsDown, secs, keyCode] = KbCheck(devNum);
        end
        if keyIsDown
            rtx = secs - timeFix; % measure from start of trial
            key = find(keyCode);
            extraKeyPresses = [extraKeyPresses; secs rtx key(1)];
                        
            if p.triggersOn
                PTBSendTrigger(p.triggers.response, 0);
                triggerTimes = [triggerTimes; p.triggers.response GetSecs];
            elseif p.triggerEEG
                ppdev_mex('Write', 1, p.triggers.response);
                WaitSecs(0.005);
                ppdev_mex('Write', 1, 0);
                triggerTimes = [triggerTimes; p.triggers.response GetSecs];
            end
        end
    end
    
    % Present fixation
    Screen('FillRect', window, white*p.backgroundColor);
    drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders, p.ringColor*white)
    drawFixation(window, cx, cy, p.fixSize, p.fixColor*white);
    if strcmp(p.testingLocation,'MEG') || strcmp(p.testingLocation,'denisonlab-EEG')
        drawPhotodiode(window, [cx cy]*2, white, 0); % black
    end
    if p.triggerEEG
        ppdev_mex('Write', 1, p.triggers.fixation); 
    end
    timeFix = Screen('Flip', window, timeIm + p.gratingDur - slack);
    if p.triggerEEG
        ppdev_mex('Write', 1, 0); 
        triggerTimes = [triggerTimes; p.triggers.fixation GetSecs];
    end
    if p.triggersOn
        PTBSendTrigger(p.triggers.fixation, 0);
        triggerTimes = [triggerTimes; p.triggers.fixation GetSecs];
    elseif p.triggerEEG
%         ppdev_mex('Write', 1, p.triggers.fixation);
%         WaitSecs(0.005);
%         ppdev_mex('Write', 1, 0);
%         triggerTimes = [triggerTimes; p.triggers.fixation GetSecs];
    end
    if p.eyeTracking
        Eyelink('Message', 'EVENT_FIX');
    end
    
    % Collect response
    responseKey = [];
    while isempty(responseKey) && GetSecs < timeIm + p.responseWindowDur
        if p.useKbQueue
            [keyIsDown, firstPress] = KbQueueCheck();
            if keyIsDown
                secs = min(firstPress(firstPress~=0));
                keyCode = firstPress==secs;
            else
                keyCode = [];
            end
        else
            [keyIsDown, secs, keyCode] = KbCheck(devNum);
        end
        responseKey = find(p.keyCodes==find(keyCode)); % must press a valid key
        if numel(responseKey)>1 % if more than one key was pressed simultaneously
            responseKey = 1;
        end
    end
    if isempty(responseKey)
        response = 0; % absent
        responseKey = NaN;
        timeResp = NaN;
        rt = NaN;
    else
        response = 1; % present
        timeResp = secs;
        rt = timeResp - timeIm; % measured from image onset
        
        if p.triggersOn
            PTBSendTrigger(p.triggers.response, 0);
            triggerTimes = [triggerTimes; p.triggers.response GetSecs];
        elseif p.triggerEEG
            ppdev_mex('Write', 1, p.triggers.response);
            WaitSecs(0.005);
            ppdev_mex('Write', 1, 0);
            triggerTimes = [triggerTimes; p.triggers.response GetSecs];
        end
        
        if p.eyeTracking
            Eyelink('Message', 'EVENT_RESPONSE');
        end
    end
    
    % Feedback
    if response==targetState
        correct = 1;
        timeTone = NaN;
        toneIdx = NaN;
    else
        correct = 0;
        if response==0 % miss
            toneIdx = find(strcmp(p.toneNames,'miss'));
            toneRefTime = timeIm + p.responseWindowDur; % end of response window
        elseif response==1 % false alarm
            toneIdx = find(strcmp(p.toneNames,'fa'));
            toneRefTime = timeResp;
        end
        feedbackTone = [p.tones(toneIdx,:); p.tones(toneIdx,:);];
        
        % Present feedback tone
        PsychPortAudio('FillBuffer', pahandle, feedbackTone*p.soundAmp);
        timeTone = PsychPortAudio('Start', pahandle, [], toneRefTime + p.toneOnsetSOA, 1);
        if p.triggersOn
            PTBSendTrigger(p.triggers.tone, 0);
            triggerTimes = [triggerTimes; p.triggers.tone GetSecs];
        elseif p.triggerEEG
            ppdev_mex('Write', 1, p.triggers.tone);
            WaitSecs(0.005);
            ppdev_mex('Write', 1, 0);
            triggerTimes = [triggerTimes; p.triggers.tone GetSecs];
        end
        if p.eyeTracking
            Eyelink('Message', 'EVENT_TONE');
        end
    end
    
    % If there was a response or a tone, wait before starting the next trial
    if response==1 || correct==0
        extraITI = p.extraITI;
    else
        extraITI = 0;
    end
    
    if p.eyeTracking
        Eyelink('Message', 'TRIAL_END');
    end
    
    % Adjust staircase level, only if target was present
    if p.staircase && targetState==1
        [stairIdx, lastFewAcc] = updateStaircase(...
            p.stairs, stairIdx, lastFewAcc, correct);
        stairValues(stairCounter) = stairIdx;
        stairCounter = stairCounter + 1;
        
        % Show the current threshold estimate
        reversalValues = getReversalValues(stairValues);
        if numel(reversalValues)>=5
            % average the last 5 reversals to get threshold
            threshold = mean(p.stairs(reversalValues(end-4:end)));
            fprintf('Threshold estimate = %f\n', threshold)
        else
            threshold = [];
        end
    end
    
    % Store trial info
    trials(iTrial,iti2Idx) = iti; % the actual ITI, including extra delay
    trials(iTrial,fixColorIdx) = fixColor(1);
    trials(iTrial,responseKeyIdx) = responseKey;
    trials(iTrial,responseIdx) = response;
    trials(iTrial,correctIdx) = correct;
    trials(iTrial,rtIdx) = rt;
    trials(iTrial,spatialLocationIdx) = spatialLocation;
    
    % Store timing
    timing.timeFix(iTrial,1) = timeFix;
    timing.timeIm(iTrial,1) = timeIm;
    timing.timeTone(iTrial,1) = timeTone;
    timing.timeResp(iTrial,1) = timeResp;
    
    % Store presented trial info
    trialsPresented(iTrial).iti = iti;
    trialsPresented(iTrial).orientation = orientation;
    trialsPresented(iTrial).phase = phase;
    trialsPresented(iTrial).contrast = contrast;
    trialsPresented(iTrial).targetState = targetState;
    trialsPresented(iTrial).fixColor = fixColor(1);
    trialsPresented(iTrial).toneIdx = toneIdx;
    trialsPresented(iTrial).extraKeyPresses = extraKeyPresses;
    trialsPresented(iTrial).spatialLocation = spatialLocation;

    if mod(iTrial,p.nTrialsPerBlock)==0 || iTrial==nTrials
        % Save the workspace every block
        save('data/TEMP')
        
        % Calculate block accuracy
        blockStartTrial = (iTrial/p.nTrialsPerBlock)*p.nTrialsPerBlock - p.nTrialsPerBlock + 1;
        if blockStartTrial < 0 % we are doing less than one block
            blockStartTrial = 1;
        end
        trialsInBlock = trials(blockStartTrial:iTrial,:);
        ts = trialsInBlock(:,targetStateIdx);
        blockHit = mean(trialsInBlock(ts==1,correctIdx));
        blockFA = 1-mean(trialsInBlock(ts==2,correctIdx));
        
        accMessage = sprintf('Hit rate: %d%%\nFalse alarm rate: %d%%', ...
            round(blockHit*100), round(blockFA*100));
        blockMessage = sprintf('%s You''ve completed %d of %d blocks.', highpraise, block, ceil(nTrials/p.nTrialsPerBlock));
        if iTrial==nTrials
            keyMessage = '';
        else
            keyMessage = 'Press 1 to go on.';
        end
        fprintf('Block %i completed. %s', block, accMessage);

        breakMessage = sprintf('%s\n%s\n\n%s', blockMessage, accMessage, keyMessage);
        DrawFormattedText(window, breakMessage, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window);
        WaitSecs(1);
        if iTrial < nTrials
            keyPressed = 0;
            while ~keyPressed
                if p.useKbQueue
                    [keyIsDown, firstPress] = KbQueueCheck();
                    keyCode = logical(firstPress);
                else
                    [secs, keyCode] = KbWait(devNum);
                end
                
                if strcmp(KbName(keyCode),'1!')
                    keyPressed = 1;
                end
            end
        end
        
        block = block+1; % keep track of block for block message only
        
        if iTrial<nTrials
            % Present initial fixation
            drawPlaceholders(window, white, p.backgroundColor*white, phRect, p.phLineWidth, p.showPlaceholders, p.ringColor*white)
            drawFixation(window, cx, cy, p.fixSize, p.fixColor*white);
            timeFix = Screen('Flip', window);
        end
    end
end
timing.endTime = GetSecs;

WaitSecs(2);
DrawFormattedText(window, 'All done! Thanks for your effort!', 'center', 'center', [1 1 1]*white);
Screen('Flip', window);
WaitSecs(2);

%% Calculate more timing things
timing.triggers = triggerTimes;
timing.timeFix = [timing.startTime; timing.timeFix];
timing.imDur = timing.timeFix(2:end) - timing.timeIm; % fix2 - im1
timing.iti = timing.timeIm - timing.timeFix(1:end-1); % im1 - fix1, trial is iti then im
timing.respToneSOA = timing.timeTone - timing.timeResp;
timing.imToneSOA = timing.timeTone - timing.timeIm;

%% Store staircase 
staircase.stairValues = stairValues;
staircase.reversalValues = reversalValues;
staircase.threshold = threshold;

%% Store expt info
expt.subjectID = subjectID;
expt.run = run;
expt.p = p;
expt.timing = timing;
expt.trials_headers = trials_headers;
expt.trials = trials;
expt.trialsPresented = trialsPresented;
expt.staircase = staircase;
expt.whenSaved = datestr(now);

%% Analyze and save data
save('data/TEMP')

if saveData
    save(dataFile, 'expt')
end

%% Save eye data and shut down the eye tracker
if p.eyeTracking
    rd_eyeLink('eyestop', window, {eyeFile, eyeDataDir});
    
    % rename eye file
    copyfile(sprintf('%s/%s.edf', eyeDataDir, eyeFile), eyeFileFull)
    delete(sprintf('%s/%s.edf', eyeDataDir, eyeFile))
end

%% Plot figures
if plotTimingFigs
    plotTiming(expt, saveTimingFigs)
end

if plotStaircaseFigs
    plotStaircase(expt, saveStaircaseFigs)
end

%% Clean up
ShowCursor;
PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
if p.triggerEEG
    ppdev_mex('Close',1);
end
Screen('CloseAll')
