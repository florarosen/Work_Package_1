% Creates audio stimuli for XPTempo experiment 

clear

savePath = fullfile('.','/stimuli');

if ~isdir(savePath)
    mkdir(savePath)
end

% we will save all the stimuli in a big structure
par.tracks = struct();

% first rhythm
par.tracks(1).name = 'unsyncopated';
par.tracks(1).pattern = [1 1 1 0 1 1 1 0 1 1 0 0];

% second rhythm
par.tracks(2).name = 'syncopated';
par.tracks(2).pattern = [1 1 1 1 0 1 1 1 0 0 1 0];

% sampling rate
par.fs           = 44100;

par.gridIOI      = 0.2; 

% time between two successive events (either sound or silence)
par.eventdur     = 0.2;

% duration of linear onset ramp for the sound event
par.rampon       = 0.010;

% duration of linear offset ramp for the sound event
par.rampoff      = 0.050;

% how many times the rhythmic pattern repeats in each trial
par.ncycles      = 25; % 17 cycles is 40.8s / 25 cycles is 60s(for 2.4s cycle)

% calculate duration of one long trial in seconds
par.trialdur    = par.ncycles * length(par.tracks(1).pattern) * par.gridIOI;

par.IS_NOISE    = false; % TRUE -> use white noise carrier | FALSE -> use tone

% carrier f0
par.f0 = 86;




%% synthesis

% make time vector for one trial
t = [0 : 1/par.fs : par.trialdur-1/par.fs];

% make carrier for one trial
if par.IS_NOISE
    carrier = rand(size(t));
else
    carrier = sin(2*pi*t*par.f0);
end

% make sure there is no clipping
carrier = carrier .* max(abs(carrier));


% make envelope of one sound event
envEvent = ones(1, round(par.eventdur * par.fs));

% apply onset and offset ramp
envEvent(1:round(par.rampon*par.fs)) = envEvent(1:round(par.rampon*par.fs)) .* linspace(0,1,round(par.rampon*par.fs));

envEvent(end-round(par.rampoff*par.fs)+1:end) = envEvent(end-round(par.rampoff*par.fs)+1:end) .* linspace(1,0,round(par.rampoff*par.fs));


% go over rhythms
for rhythmi=1:length(par.tracks)

    % allocate envelope vector for the whole trial (as zeros)
    envTrial = zeros(size(carrier));

    % make ncycles copies of the rhythmic pattern
    patternTrial = repmat(par.tracks(rhythmi).pattern, 1, par.ncycles);

    % go over each event in the trial
    for i=1:length(patternTrial)
        % if the event is sound
        if patternTrial(i)
            % find the time of event onset
            eventTime = (i-1)*par.gridIOI;
            % convert to index
            eventIdx = round(eventTime*par.fs);
            % paste the sound event envelope
            envTrial(eventIdx+1:eventIdx+length(envEvent)) = envEvent;
        end
    end

    % multiply carrier and envelope for the whole trial
    s = carrier .* envTrial;
    %s = s';

    % make it stereo
    s = [s',s'];

    % save into structure
    par.tracks(rhythmi).env = envTrial;
    par.tracks(rhythmi).s = s;

    % also write wav file
    audiowrite(fullfile(savePath,sprintf('%ds_%dhz_%s.wav',par.trialdur,par.f0,par.tracks(rhythmi).name)), s, par.fs);
end

% save the structure
save(fullfile(savePath,sprintf('%ds_%dhz_stimuli.mat',par.trialdur,par.f0)), 'par');