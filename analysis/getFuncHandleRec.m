% Look up the handle for a function which works on a single recording.
%
% The function handle is taken from a subfunction within this .m file.
% All functions returned have the format
%     [X,W] = func(signal, tt, auxlab, ped, vel, acc)
% for stimulus-independent or
%     [X,W] = func(signal, tt, auxlab, ped, vel, acc, stimTimes)
% for stimulus dependent functions.
%
% All inputs should be arrays sized [1, numTimePoints], except
% `signal`, which should be [numROIs, numTimePoints]. Here numROIs is
% the number of non-trivial ROIs present in the recordings, numTimePoints
% is the number of samples over the duration of any single recording (the
% duration in seconds times the sampling frequency).
% The output X is sized [numROIs, ...], with the size of higher dimensions
% specific to the function called, and contains the output of the function.
% W is sized the same as X and, unless otherwise specified, contains the count 
% or duration over which X was taken for this recording. If all X should 
% have the same weighting, this will be an array of ones.
%
% Inputs
% ------
% funcstr : string
%     Name of the function.
%
% Outputs
% -------
% func : function handle
%     Handle to the requested function.
%
% ======================================================================= %
%                           LIST OF FUNCTIONS                             %
% ======================================================================= %
% (Ctrl+D to jump to function in this file)
%__________________________________________________________________________
%   # Functions using VRdata
%   averageVRbyTrackRwds
%   averageVRbyTrackLicksRwdNonrwd
%   TestLicksVSnonLickPreRwd
%   VRsuccessRate
%__________________________________________________________________________


% ======================================================================= %
%                           MAIN FUNCTION                                 %
% ======================================================================= %

function func = getFuncHandleRec(funcstr)

% Look up the function by name from the local scope
% We find the function with the correct name from below and return a handle
% to it!
% But take caution here; str2func always succeeds, even if the function
% doesn't exist.
func = str2func(funcstr);

% Make sure the function actually exists
try
    func(); % This is expected to break due to lack of input
    % If it did not break, raise a warning 
    warning('ROCHLAB:functionWithoutInputs', ...
        'Function %s ran without any inputs', funcstr);
catch ME
    if strcmp(ME.identifier,'MATLAB:UndefinedFunction') || ...
        strcmp(ME.message, ['Undefined function or variable ''' funcstr ''''])
        % If we cannot call the function, rethrow the error
        rethrow(ME);
    end
    % If it was any other error, we don't mind
end
    
end


% ======================================================================= %
%                       FUNCTIONS USING VRdata                            %
% ======================================================================= %
% Description of possible inputs
%__________________________________________________________________________
% Inputs
% ------
%   signal : array [nROIs, nSamples]
%       Fluorescence activity traces per ROI for each trial and track
%       nSamples = tot number of samples during the recording
%   tt : vector [1, nSamples]
%       Time of samples since the start of the recording (seconds)
%   auxlab -> here, VRlab : structure with fields:
%       Each field is a vector of size [1, nSamples]
%       - trackNum (index of VR track presented by time/frame)
%       - reward (event detection of reward; 0= no reward, 1 = early reward, 2 = default (late) reward)
%       - trialNum (running index of trial number by time/frame)
%       - licks (even detection of licking from lick sensor; 0 = no lick, 1 = lick)
%       - actlab (index of 'action' of animal by time/frame; 0 = stationary; 3 = locomoting)
%__________________________________________________________________________


function [X, W] = averageVRbyTrackRwds(signal, tt, auxlab, ~, ~, ~, ~)
% For each track, extracts the average DF/F0 signal of each ROI for various
% reward types (labels of rewards are given below).
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nRewardTypes+1, 1, nTracks] 
%           Average signal.
%   W : array [nROIs, nRewardTypes+1, 1, nTracks]
%           For each track, returns the number of samples of each reward.
%           Note that all rows except the first are just NaN or copies.
%   NOTE: The last data (nRewardTypes+1) consider all reward types together
%__________________________________________________________________________


[dataDf_tracks, dataReward_tracks] = prepdataforfunction_Rec('averageVRbyTrackRwds', signal, tt, auxlab);
% Cell arrays [1, nTracks], where for each track
%   dataDf_tracks : the DF/F0 of each ROI [nROIs, nSamples_track]
%   dataReward_tracks : the labels of reward type [1, nSamples_track]

rewardTypes_list = 0:4;
% List of rew types:
%       0 = track(nonRewardtime)
%       1 = reward
%       2 = defaultRewd
%       3 = preReward
%       4 = preDefaultRwd

nTracks = length(dataDf_tracks);
nRewardTypes = length(rewardTypes_list);
nROIs = size(signal,1);

X = nan(nROIs, nRewardTypes+1, 1, nTracks);
W = nan(nROIs, nRewardTypes+1, 1, nTracks);
% Take average from each action type
for iTrack = 1:nTracks
    if ~isempty(dataDf_tracks{iTrack})
        
        for iAuxlab = 1:nRewardTypes
            is_sample_reward = dataReward_tracks{iTrack} == rewardTypes_list(iAuxlab);
            X(:,iAuxlab,1,iTrack) = mean(dataDf_tracks{iTrack}(:,is_sample_reward),2);
            W(1,iAuxlab,1,iTrack) = sum(double(is_sample_reward));
        end
        X(:,nRewardTypes+1,1,iTrack) = mean(dataDf_tracks{iTrack},2);
        W(1,nRewardTypes+1,1,iTrack) = length(dataReward_tracks{iTrack});
    end
end

end

function [X, W] = averageVRbyTrackLicksRwdNonrwd(signal, tt, auxlab, ~, ~, ~, ~)
% For each track, extracts the average DF/F0 signal of each ROI for various
% action types (labels of actions are given below).
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nActionTypes, 1, nTracks] (average signal)
%   W : array [nROIs, nActionTypes, 1, nTracks]
%           For each track, returns the number of samples of the action.
%           Note that all rows are just copies.
%__________________________________________________________________________

avgHz = 1/(mean(diff(tt)));
frames = floor(avgHz*2); % window around (before and after) lick that gets counted as lick (since it is normally a single frame event


[dataDf_tracks, dataReward_tracks, dataLicks_tracks] = prepdataforfunction_Rec('averageVRbyTrackLicksRwdNonrwd', signal, tt, auxlab);
% Cell arrays [1, nTracks], where for each track
%   dataDf_tracks : the DF/F0 of each ROI [nROIs, nSamples_track]
%   dataReward_tracks : the labels of reward type [1, nSamples_track]
%   dataLicks_tracks : the labels lick / no lick [1, nSamples_track]

pairs_LicksRewardLabel = struct(...
    'licksLabel',   {0, 1, 1, 1, 0, 0}, ...
    'rewardLabel',  {0, 0, 1, 2, 1, 2});
% List of labels for actions:
%       Licks   Reward
%   1)  0       0       no lick pre-reward
%   2)  1       0       licks pre-reward
%   3)  1       1       licks post hit reward
%   4)  1       2       licks post miss reward
%   5)  0       1       no lick post hit reward
%   6)  0       2       no lick post miss reward

nTracks = length(dataDf_tracks);
nActionTypes = length(pairs_LicksRewardLabel);
nROIs = size(signal,1);
X = nan(nROIs, nActionTypes, 1, nTracks);
W = nan(nROIs, nActionTypes, 1, nTracks);

% Take average from each action type
for iTrack = 1:nTracks
    
    if ~isempty(dataDf_tracks{iTrack})
        
        for iAuxlab = [2,3,4,6]
            is_sample_action = ( dataReward_tracks{iTrack} == pairs_LicksRewardLabel(iAuxlab).rewardLabel ) & ...
                ( dataLicks_tracks{iTrack} == pairs_LicksRewardLabel(iAuxlab).licksLabel );
            
            X(:,iAuxlab,1,iTrack) = mean(dataDf_tracks{iTrack}(:,is_sample_action),2);
            W(:,iAuxlab,1,iTrack) = sum(double(is_sample_action));
        end
        
        for iAuxlab = [1,5]
            is_sample_action = ( dataReward_tracks{iTrack} == pairs_LicksRewardLabel(iAuxlab).rewardLabel ) & ...
                ( dataLicks_tracks{iTrack} == pairs_LicksRewardLabel(iAuxlab).licksLabel );
            
            if iAuxlab==1
                LickCnt = W(1,2,1,iTrack);
            else
                LickCnt = W(1,3,1,iTrack);
            end
            
            if sum(double(is_sample_action)) > LickCnt*2
                for iBlock = 1:length(dataReward_tracks{iTrack})/(frames*2)
                    is_sample_action(1,iBlock*frames*2:(iBlock*frames*2)+frames) = false;
                end
            end
            X(:,iAuxlab,1,iTrack) = mean(dataDf_tracks{iTrack}(:,is_sample_action),2);
            W(:,iAuxlab,1,iTrack) = sum(double(is_sample_action));
            
        end 
    end
end

end

function [X, W] = TestLicksVSnonLickPreRwd(signal, ~, auxlab, ~, ~, ~, ~)
% For each track, for each ROI, compare signal when animal licks versus not
% before reward. Nonparametric Kruskal–Wallis one-way analysis of variance.
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nTracks] (Kruskal–Wallis test)
%   W : logical array [nROIs, nTracks]
%           For each ROI, returns whether average signal when licks is
%           higher than when does not lick.
%__________________________________________________________________________

[dataDf_tracks, dataReward_tracks, dataLicks_tracks] = prepdataforfunction_Rec('TestLicksVSnonLickPreRwd', signal, auxlab);
% Cell arrays [1, nTracks], where for each track
%   dataDf_tracks : the DF/F0 of each ROI [nROIs, nSamples_track]
%   dataReward_tracks : the labels of reward type [1, nSamples_track]
%   dataLicks_tracks : the labels lick / no lick [1, nSamples_track]

alpha = 0.001; % significance level for p-value

nROIs = size(signal,1);
nTracks = length(dataDf_tracks);

X = nan(nROIs,nTracks); % see below for column designations
W = nan(nROIs,nTracks);

for iTrack = 1:nTracks
       
    if isempty(dataDf_tracks{iTrack})
        continue
    end
    
    lick_noreward = (dataLicks_tracks{iTrack}==1) & (dataReward_tracks{iTrack}==0);
    lickDataDf = dataDf_tracks{iTrack}(:,lick_noreward);
    nLicks = size(lickDataDf,2);
    
    nolick_noreward = (dataLicks_tracks{iTrack}==0) & (dataReward_tracks{iTrack}==0);
    % Limit non-lick data to same amount as lick data
    spfactor = floor(sum(nolick_noreward)/nLicks);
    if  sum(nolick_noreward)>nLicks*spfactor
        for iBlock = 1:size(nolick_noreward,2)/spfactor
            nolick_noreward(randperm(size(lnolick_noreward,2)-floor(nLicks*1.1)))=0;
        end
    end
    nonLickDataDf = dataDf_tracks{iTrack}(:,nolick_noreward);
   
    for iROI = 1:nROIs
        X(iROI,iTrack) = signrank(lickDataDf(iROI,:), nonLickDataDf(iROI,:), 'alpha', alpha);
        W(iROI,iTrack) = nanmean(lickDataDf(iROI,:)) > nanmean(nonLickDataDf(iROI,:));
    end
end


end

function [X, W] = VRsuccessRate(signal, ~, auxlab, ~, ~, ~, ~)
% Returns success rate for each track (successful trial = hit reward)
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 3, 1, nTracks+1]
%       columns :   1 - success rate
%                   2 - number of successful trials
%                   3 - number of unsuccessful trials
%   W : array [nROIs, 3, 1, nTracks+1]
%       columns :   1 - number of successful trials
%                   2 - number of unsuccessful trials
%                   3 - number of trials
%   Note that all rows in X and W are just copies as sucess rate applies to all ROIs.
%__________________________________________________________________________


[numTrials_tracks, numhitTrials_tracks] = prepdataforfunction_Rec('VRsuccessRate', auxlab);
% numTrials_tracks : array [1, nTracks]
% numhitTrials_tracks : array [1, nTracks]
%   NOTE: if there is no data for the track, the entry is NaN

nROIs = size(signal,1);
nTracks = length(numTrials_tracks);

X = nan(1, 3, 1, nTracks+1);
W = nan(1, 3, 1, nTracks+1);

% Percent success rate
X(1,1,1,1:nTracks) = numhitTrials_tracks ./ numTrials_tracks; 
% Num hits
X(1,2,1,1:nTracks) = numhitTrials_tracks;
W(1,1,1,1:nTracks) = numhitTrials_tracks; 
% Num misses
X(1,3,1,1:nTracks) = numTrials_tracks - numhitTrials_tracks;
W(1,2,1,1:nTracks) = X(1,3,1,1:nTracks); 
% Number of trials
W(1,3,1,1:nTracks) = numTrials_tracks; 

% Same, but for all tracks together
numTrials = sum(numTrials_tracks);
numhitTrials = sum(numhitTrials_tracks);
% Overall success rate
X(1,1,1,end) = numhitTrials / numTrials; 
% Num hits
X(1,2,1,end) = numhitTrials;
W(1,1,1,end) = numhitTrials; 
% Num misses
X(1,3,1,end) = numTrials - numhitTrials;
W(1,2,1,end) = X(1,3,1,end); 
% Number of trials
W(1,3,1,end) = numTrials; 

X = repmat(X, [nROIs, 1, 1, 1]);
W = repmat(W, [nROIs, 1, 1, 1]);

end