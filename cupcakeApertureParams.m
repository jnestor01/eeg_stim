function p = cupcakeApertureParams

p.testingLocation = 'denisonlab-EEG'; %'desk'; % 'CarrascoL1',, 'MEG'
KbName('UnifyKeyNames');
switch p.testingLocation
%     case 'MEG'
%         p.keyNames = {'1!'};
%         p.refRate = 60;
%         p.screenSize = [23.5 17]; % cm
%         p.screenRes = [1024 768];
%         p.viewDist = 42; % cm
%         p.eyeTracking = 1;
%         p.useKbQueue = 1;
%         p.soundAmp = 0.1;
%         p.triggersOn = 1;
%         p.displayPath = '/Users/megadmin/Desktop/Experiments/Rachel/vistadisp/exptTools2/displays/meg_lcd_20180420_brightness-32';
    case {'desk'}
        p.keyNames = {'1!'};
        p.refRate = 60;
        p.screenSize = [13 9]; % (in)
        p.screenRes = [1280 1024];
        p.viewDist = 36; % (in)
        p.eyeTracking = 0;
        p.useKbQueue = 0;
        p.soundAmp = 1;
        p.triggersOn = 0;
%     case 'CarrascoL1'
%         p.keyNames = {'1!'};
%         p.refRate = 60;
%         p.screenSize = [40 30];
%         p.screenRes = [1280 960];
%         p.viewDist = 56;
%         p.eyeTracking = 0; 
%         p.useKbQueue = 0;
%         p.soundAmp = 1;
%         p.triggersOn = 0;
    case {'denisonlab-EEG'}
        p.keyNames = {'1!'};
        p.refRate = 120;
        p.screenSize = [53.5 30]; % (in)
        p.screenRes = [1280 1024];
        p.viewDist = 57; % (in)
        p.eyeTracking = 0;
        p.useKbQueue = 0;
        p.soundAmp = 1;
        p.triggersOn = 0;
        p.triggerEEG = 1;
    otherwise
        error('Testing location not found.')
end

p.keyCodes = KbName(p.keyNames);
p.backgroundColor = 0.35;
p.nReps = 1; % 6 reps with 32 trials/rep -> 3 blocks
p.nTrialsPerBlock = 60; %64 trials -> ~1.5 min/block
p.eyeRad = 4; % allowed fixation radius (degrees)   

% Text
p.font = 'Verdana';
p.fontSize = 24;

% Placeholders
% 2 for circle placeholders
p.showPlaceholders = 2;
p.phLineWidth = 4; % (pixels)
p.phCushion = 0.75; % degrees, expand ph by this amount relative to the image

p.ringColor = 0.5; % color inside region framed by ring placeholders

% Fixation
p.fixColor = 1; %[1; p.backgroundColor]; % [inner; outer] 
% p.fixDiameter = .35; %[.35 .7]; % deg
p.fixDiameter = .85; %[.35 .7]; % deg

% Timing
p.itiType = 'hazard'; % 'uniform','hazard'
switch p.itiType
    case 'uniform'
        p.itis = 1:0.05:1.8;
    case 'hazard'
        p.itis = 1:0.05:2.5;
        p.hazardProb = 0.2; % at every time step, the probability of the event is 0.2
        p.itiPDF = p.hazardProb.*(1-p.hazardProb).^(0:numel(p.itis)-1); % f = p.*(1-p).^x;
    otherwise
        error('p.itiType not recognized')
end
p.extraITI = 0.5; % insert after a button press or feedback tone

p.gratingDur = 0.5; 
p.toneDur = 0.2;
p.responseWindowDur = 0.5;
p.eyeSlack = 0.5; % minimum time between regain fixation and next stimulus

% Images
p.imPos = [0 0];
% p.imSize = [17 17]; % this is the size of the image container that holds the stim
% p.gratingDiameter = [16 1]; % [outer inner] 
p.imSize = 3.5;
p.gratingDiameter = [3 0]; % [outer inner] 

p.gratingSF = 1.5; % cpd
% p.gratingOrientations = [0:179]; %0:22.5:179; %0:5:179; %0:20:179; %0:22.5:179
% p.gratingOrientations = [0:5:175]; %0:22.5:179; %0:5:179; %0:20:179; %0:22.5:179
% p.gratingOrientations = [0:45:135];
p.gratingOrientations = 0;

% number of evenly distributed spatial locations, starting at vertical
% set radial distance to 0 for simple orientation task
p.spatialLocations = linspace(0,360,36+1);
p.spatialLocations(end) = [];
p.radialDistance = 5.5;

p.gratingPhases = 0; % eg. 0, or [0 pi/2 pi 3*pi/2]
p.gratingContrasts = [0.25 1]; 
p.noiseContrasts = 0; 
p.aperture = 'cosine'; % disk:'cosine', annulus:'cosine-ring', concentric:'radial-sine-ring', roth:'vignette-ring'
p.apertureEdgeWidth = 1; % half of a period, so sf of radial-sine aperture is 1/(2*width)
p.angularFreq = 8; % for 'vignette' only
if strfind(p.aperture,'radial')
    p.apertureSF = 1/(2*p.apertureEdgeWidth);
end

% Task
p.targetStates = [1 0]; % 1=present, 0=absent
p.propTargetPresent = .2;
% add orientation 

% Staircase
p.staircase = 1;
p.stairs = 1-logspace(-1.5,0,15); % hard (~1) to easy (0)
if p.staircase
    p.targetContrast = 0;
    fprintf('\nStaircase is ON\n')
end

% Sounds
p.Fs = 44100;
p.toneNames = {'miss','fa'};
p.toneFreqs = [523 784]; % [lower C, higher G]
for iTone = 1:numel(p.toneFreqs)
    tone = MakeBeep(p.toneFreqs(iTone), p.toneDur, p.Fs);
    p.tones(iTone,:) = applyEnvelope(tone, p.Fs);
end
p.toneOnsetSOA = 0.02; % 20 ms
% 10^0.5 for every 10dB 

% MEG triggers
p.triggers.fixation = 2^0; %1
p.triggers.image = 2^1; %2
p.triggers.tone = 2^2; %4
p.triggers.target = 2^3; %8
p.triggers.response = 2^4; %16

% set a seed to replay a past run
p.setSeed = [];
rng("shuffle");
p.randSeed = rng().Seed;

