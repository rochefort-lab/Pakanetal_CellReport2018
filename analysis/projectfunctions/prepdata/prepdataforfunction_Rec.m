function [varargout] = prepdataforfunction_Rec(nameOfFunction, varargin)


% *********************************************************************
% -- VARIOUS FUNCTIONS USING VR data
% *********************************************************************
switch nameOfFunction
    
    case 'averageVRbyTrackRwds'
        
        signal_df = varargin{1};
        tt = varargin{2};
        VRlab = varargin{3};
        
        [dataDf_tracks, dataReward_tracks] = prepare_averageVRbyTrackRwds(signal_df, tt, VRlab);
        
        varargout{1} = dataDf_tracks;
        varargout{2} = dataReward_tracks;
        
        
    case 'averageVRbyTrackLicksRwdNonrwd'
        
        signal_df = varargin{1};
        tt = varargin{2};
        VRlab = varargin{3};
        
        [dataDf_tracks, dataReward_tracks, dataLicks_tracks] = prepare_averageVRbyTrackLicksRwdNonrwd(signal_df, tt, VRlab);
        
        varargout{1} = dataDf_tracks;
        varargout{2} = dataReward_tracks;
        varargout{3} = dataLicks_tracks;
        
        
    case 'pairedTestRwdsBinbyTime'
        
        signal_df = varargin{1};
        tt = varargin{2};
        VRlab = varargin{3};
        numTracks = 8;
        
        [dataDf_PostReward, dataDf_PreReward] = prepare_pairedTestRwdsBinbyTime(signal_df, tt, VRlab);
        
        varargout{1} = dataDf_PostReward;
        varargout{2} = dataDf_PreReward;
        varargout{3} = numTracks;
        
        
    case 'TestLicksVSnonLickPreRwd'
        
        signal_df = varargin{1};
        VRlab = varargin{2};
        
        [dataDf_tracks, dataReward_tracks, dataLicks_tracks] = prepare_TestLicksVSnonLickPreRwd(signal_df, VRlab);
        
        varargout{1} = dataDf_tracks;
        varargout{2} = dataReward_tracks;
        varargout{3} = dataLicks_tracks;
        
        
    case 'VRsuccessRate'
        
        VRlab = varargin{1};
        
        [numTrials_tracks, numhitTrials_tracks] = prepare_VRsuccessRate(VRlab);
        
        varargout{1} = numTrials_tracks;
        varargout{2} = numhitTrials_tracks;
        
end




end

%__________________________________________________________________________
% ___ FUNCTIONS ___________________________________________________________


function [dataDf_tracks, dataReward_tracks] = prepare_averageVRbyTrackRwds(signal_df, tt, VRlab)

avgHz = 1/(mean(diff(tt)));
frames = floor(avgHz);
dataReward = VRlab.reward;
numSamples = length(dataReward);

% First expand actions to 1000ms after signal (since rewards and licks are
% single frame event recordings)
% NOTE: 0 = track(nonRewardtime); 1 = reward; 2 = defaultRewd
list_indices = 1:numSamples;
is_possible_startIndex = list_indices >= frames+1;
for iRewardType = 1:2
    % Find reward of this type and make sure that an event was not too soon
    is_sample_reward = (VRlab.reward == iRewardType) & is_possible_startIndex;
    
    if any(is_sample_reward)
        Idx_reward = find(is_sample_reward);
        for iEvent = 1:length(Idx_reward)
            dataReward(1,Idx_reward(iEvent):Idx_reward(iEvent)+frames) = iRewardType;
            dataReward(1,Idx_reward(iEvent)-(frames+1):Idx_reward(iEvent)) = iRewardType+2;
        end
    end
end
% Make sure we did not add data points on end
dataReward = dataReward(1,1:size(tt,2));

% Organise per track
numTracks = 8;
dataReward_tracks = cell(1,numTracks);
dataDf_tracks = cell(1,numTracks);
for iTrack = 1:numTracks
    is_sample_track = VRlab.trackNum==iTrack;
    dataReward_tracks{iTrack} = dataReward(is_sample_track);
    dataDf_tracks{iTrack} = signal_df(:, is_sample_track);
end


end

function [dataDf_tracks, dataReward_tracks, dataLicks_tracks] = prepare_averageVRbyTrackLicksRwdNonrwd(signal_df, tt, VRlab)

avgHz = 1/(mean(diff(tt)));
frames = floor(avgHz*2); % window around (before and after) lick that gets counted as lick (since it is normally a single frame event

dataReward = VRlab.reward;
dataLicks = VRlab.licks;
numSamples = length(dataLicks);
nTrials = max(VRlab.trialNum);

% First expand actions to 500ms before and 500ms after signal (since
% rewards and licks are single frame event recordings)
iLicksLab = 1;
list_indices = 1:numSamples;
is_possible_startIndex = list_indices >= frames+1;
% Find data of this type and make sure that an event was not too soon
is_sample_licks = (VRlab.licks == iLicksLab) & is_possible_startIndex;
if any(is_sample_licks)
    Idx_licks = find(is_sample_licks);
    for iEvent = 1:length(Idx_licks)
        dataLicks(1,Idx_licks(iEvent)-frames:Idx_licks(iEvent)+frames) = iLicksLab;
    end
end
% Make sure we did not add data points on end
dataLicks = dataLicks(1,1:size(tt,2));

for iRewardType = 1:2
    % Find reward of this type and make sure that an event was not too soon
    is_sample_reward = (VRlab.reward == iRewardType);
    
    if any(is_sample_reward)
        for iTrial = 1:nTrials
            is_sample_trial = VRlab.trialNum==iTrial;
            is_sample_reward_trial = is_sample_reward & is_sample_trial;
            if any(is_sample_reward_trial)
                Idx_start = find(is_sample_reward_trial, 1, 'first');
                Idx_stop = find(is_sample_trial, 1, 'last');
                dataReward(Idx_start:Idx_stop) = iRewardType;
            end
        end
    end
    
end

% Organise per track
numTracks = 8;
dataLicks_tracks = cell(1,numTracks);
dataReward_tracks = cell(1,numTracks);
dataDf_tracks = cell(1,numTracks);
for iTrack = 1:numTracks
    is_sample_track = VRlab.trackNum==iTrack;
    if any(is_sample_track)
        dataLicks_tracks{iTrack} = dataLicks(is_sample_track);
        dataReward_tracks{iTrack} = dataReward(is_sample_track);
        dataDf_tracks{iTrack} = signal_df(:, is_sample_track);
    end
end


end

function [dataDf_tracks, dataReward_tracks, dataLicks_tracks] = prepare_TestLicksVSnonLickPreRwd(signal_df, VRlab)

nTrials = max(VRlab.trialNum);
is_sample_reward = VRlab.reward >0;
dataReward = VRlab.reward;

for iTrial = 1:nTrials
    
    is_sample_trial = VRlab.trialNum==iTrial;
    is_sample_reward_trial = is_sample_reward & is_sample_trial;
    
    if any(is_sample_reward_trial)
        Idx_start = find(is_sample_reward_trial, 1, 'first');
        Idx_stop = find(is_sample_trial, 1, 'last');
        dataReward(Idx_start:Idx_stop) = 1;
    end
    
end

% Organise per track
numTracks = 8;
dataLicks_tracks = cell(1,numTracks);
dataReward_tracks = cell(1,numTracks);
dataDf_tracks = cell(1,numTracks);
for iTrack = 1:numTracks
    is_sample_track = VRlab.trackNum==iTrack;
    if any(is_sample_track)
        dataLicks_tracks{iTrack} = VRlab.licks(is_sample_track);
        dataReward_tracks{iTrack} = dataReward(is_sample_track);
        dataDf_tracks{iTrack} = signal_df(:, is_sample_track);
    end
end


end

function [numTrials_tracks, numhitTrials_tracks] = prepare_VRsuccessRate(VRlab)

numTracks = 8;
            
numhitTrials_tracks = nan(1, numTracks);
numTrials_tracks = nan(1, numTracks);
            

for iTrack = 1:numTracks
    
    is_sample_of_track = VRlab.trackNum == iTrack;
    
    if any(is_sample_of_track)
        
    list_trialNumbers_track = unique(VRlab.trialNum(is_sample_of_track));
    numTrials = length(list_trialNumbers_track);
    
    % Count trials that were successful (reward hit)
    numhitTrials = 0;
    for iTrial = 1:numTrials
        selected_trial = list_trialNumbers_track(iTrial);
        trial_track_samples = is_sample_of_track & (VRlab.trialNum == selected_trial);
        
        if any(VRlab.reward(trial_track_samples)==1) % Reward hit
            numhitTrials = numhitTrials+1;
        end
    end
    
    numhitTrials_tracks(iTrack) = numhitTrials; 
    numTrials_tracks(iTrack) = numTrials;
    
    end
end

end
