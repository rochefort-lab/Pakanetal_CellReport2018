% Look up the handle for a function which works on a set of recordings.
%
% The function handle is taken from a subfunction within this .m file.
% All functions returned have the format
%     X = func(signal, tt, auxlab, ped, vel, lagTime)
% for stimulus-independent or
%     [X, ustm] = func(signal, tt, auxlab, ped, vel, lagTime, stimTimes, ustm)
% for stimulus dependent functions. In the later case, the size of X in
% dimension 2 should match the number of elements in ustm, since this is
% the dimension assumed to denote different stimuli. Alternatively the
% function can return ustm as an empty array if X is not organized by
% stimulus identity.
%
% All inputs should be arrays sized [1, numTimePoints, numRecs], except
% `signal` which should be [numROIs, numTimePoints, numRecs] and both
% `stimTimes` and `ustm`, which will be cell arrays sized [numRecs, 1].
% Here numROIs is the number of non-trivial ROIs present in the recordings,
% numTimePoints is the number of samples over the duration of any single
% recording (the duration in seconds times the sampling frequency), and
% numRecs is the number of recordings to analyse together in the set. This
% input schema means the recordings must all be from the same field of view
% (so they have the same ROIs present) and must have the same duration and
% sampling frequency.
% The output X is sized [numROIs, ...], with the size of higher dimensions
% specific to the function called.
% If there is no function named funcstr within getFuncHandleRGroup.m, but
% there is one in getFuncHandleRec.m, the returned function handle is
% to an anonymous function which applies the single-recording function to
% every recording in the input and concatenates the results along the 3rd
% dimension.
%
% Inputs
% ------
% funcstr : string
%     Name of the function to lookup.
%
% Outputs
% -------
% func : function handle
%     Handle to the requested function.
%
% See also getFuncHandleRec.
%
% ======================================================================= %
%                           LIST OF FUNCTIONS                             %
% ======================================================================= %
% (Ctrl+D to jump to function in this file)
%__________________________________________________________________________
%   # Functions using VRdata
%   binnedVRData
%   binnedVRDataHit
%   binnedVRDataMiss
%   gainModVRData
%   lmiVR
%   rwdOnsetData
%   RZdistOnsetData
%   successRateSMI_VR
%   pairedTestGratVSrwdZone
%   tempMatchVR_HitMiss
%   tempMatchVR_RZ_CuedUncued
%   trialTime
%   slowMedfastTrials
%   peakVarByOnset
%   mancovaVR
%   mancovaVRshuffDist

%__________________________________________________________________________


% ======================================================================= %
%                           MAIN FUNCTION                                 %
% ======================================================================= %

function func = getFuncHandleRGroup(funcstr)

% Look up the function by name from the local scope
% We find the function with the correct name from below and return a handle
% to it!
% But take caution here; str2func always succeeds, even if the function
% doesn't exist.
func = str2func(funcstr);
try
    func(); % This is expected to break due to lack of input
    % If it did not break, raise a warning
    warning('ROCHLAB:functionWithoutInputs', ...
        'Function %s ran without any inputs', funcstr);
catch ME
    if ~strcmp(ME.identifier,'MATLAB:UndefinedFunction') && ...
            ~strcmp(ME.message, ['Undefined function or variable ''' funcstr ''''])
        % If we errored for anything other than a missing function,
        % the function is defined and we have errored due to lack of inputs
        % This is safe to return
        disp('Got function from rec group function handles');
        return;
    end
    % If the function does not exist, try to take it from the single rec
    % functions and stack the outputs together
    try
        innerfunc = getFuncHandleRec(funcstr);
        func = @(varargin) single2group(innerfunc, varargin{:});
        disp('Got function from single rec function handles, and put it in a wrapper');
        warning('ROCHLAB:unnecessarySingleFuncWrap', ...
            'Please consider using calcRecFuncPool with cat instead');
    catch ME2
        if strcmp(ME2.identifier,'MATLAB:UndefinedFunction') && ...
            strcmp(ME2.message, ['Undefined function or variable ''' funcstr '''']);
            rethrow(ME);
        else
            rethrow(ME2);
        end
    end
end

end
function out = single2group(innerfunc, varargin)

% Define dimensions
dROI  = 1;
dTime = 2;
dRec  = 3;

% Check which input args are one-for-every-rec
nArg = length(varargin);
nRecPerArg = cellfun('size', varargin, dRec);
nRec = max(nRecPerArg);

% Loop over every rec
for iRec=1:nRec
    % Assemble the subset of data to put into the single rec function
    args = cell(1,nArg);
    for iArg=1:nArg
        if nRecPerArg(iArg)==1
            args{iArg} = varargin{iArg};
        else
            args{iArg} = varargin{iArg}(:,:,iRec);
        end
    end
    % Apply the inner function to the rec
    X = innerfunc(args{:});
    % Should add some handling for X being a cell array of many outputs...
    % If this is the first rec, we need to initialise the output
    if iRec==1
        if ndims(X)>3
            error('Output has too many dimensions');
        end
        out = zeros(size(X,1),size(X,2),0);
    end
    % Stack results from each rec together
    out = cat(dRec, out, X);
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
%   disp : vector [1, nSamples]
%       Cumulative total distance travelled (displacement) since the start 
%       of the recording (cm)
%   vel : vector [1, nSamples]
%       Instantaneous velocity (cm/s)
%   lagTime : vector [1, nSamples]
%       Instantaneous frame by frame lag in VR (seconds)
%__________________________________________________________________________


function [X,stimFlag] = binnedVRData(signal, tt, auxlab, disp, vel, lagTime, ~, ~) %#ok<*DEFNU>
% Returns mean deltaf averaged across trials and binned according to track
% length (Here: 2.5cm per bin)
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, numBin_max, nTracks]
%       where numBin_max is the number of bins of the longest track
%__________________________________________________________________________

stimFlag = [];

[ meandF , numBin_max ] = prepdataforfunction_RGroup('binnedVRData',signal, tt, auxlab, disp, vel, lagTime);
%  meandF is a cell array [numTracks, 1]
%       Each cell contains an array [nROIs, numBins] where numBins is the
%       number of bins for the track.
%       The values DF/F0 were obtained by averaging over trials.
%       /!\ Tracks can be of different lengths, but here bins are fixed
%       (2.5cm) and hence numBins can be different between cells

% Number of tracks and ROIs
nTracks = length(meandF);
nROIs = size(signal,1);

% Initialise output array
X = nan(nROIs, numBin_max, nTracks); 

for iTrack = 1:nTracks
    
    if ~isempty(meandF{iTrack})
        numBins = size(meandF{iTrack}, 2);
        X(:, 1:numBins, iTrack) = meandF{iTrack};
    end
    
end

end

function [X,stimFlag] = binnedVRDataHit(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% Returns mean deltaf averaged across trials and binned according to track
% length (Here: 2.5cm per bin)
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, numBin_max, nTracks]
%       where numBin_max is the number of bins of the longest track
%__________________________________________________________________________

stimFlag = [];

[ meandF , numBin_max ] = prepdataforfunction_RGroup('binnedVRData',signal, tt, auxlab, disp, vel, lagTime,'useHitsOnly');
%  meandF is a cell array [numTracks, 1]
%       Each cell contains an array [nROIs, numBins] where numBins is the
%       number of bins for the track.
%       The values DF/F0 were obtained by averaging over trials.
%       /!\ Tracks can be of different lengths, but here bins are fixed
%       (2.5cm) and hence numBins can be different between cells

% Number of tracks and ROIs
nTracks = length(meandF);
nROIs = size(signal,1);

% Initialise output array
X = nan(nROIs, numBin_max, nTracks); 

for iTrack = 1:nTracks
    
    if ~isempty(meandF{iTrack})
        numBins = size(meandF{iTrack}, 2);
        X(:, 1:numBins, iTrack) = meandF{iTrack};
    end
    
end

end

function [X,stimFlag] = binnedVRDataMiss(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% Returns mean deltaf averaged across trials and binned according to track
% length (Here: 2.5cm per bin)
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, numBin_max, nTracks]
%       where numBin_max is the number of bins of the longest track
%__________________________________________________________________________

stimFlag = [];

[ meandF , numBin_max ] = prepdataforfunction_RGroup('binnedVRData',signal, tt, auxlab, disp, vel, lagTime,'useMissOnly');
%  meandF is a cell array [numTracks, 1]
%       Each cell contains an array [nROIs, numBins] where numBins is the
%       number of bins for the track.
%       The values DF/F0 were obtained by averaging over trials.
%       /!\ Tracks can be of different lengths, but here bins are fixed
%       (2.5cm) and hence numBins can be different between cells

% Number of tracks and ROIs
nTracks = length(meandF);
nROIs = size(signal,1);

% Initialise output array
X = nan(nROIs, numBin_max, nTracks); 

for iTrack = 1:nTracks
    
    if ~isempty(meandF{iTrack})
        numBins = size(meandF{iTrack}, 2);
        X(:, 1:numBins, iTrack) = meandF{iTrack};
    end
    
end

end

function [X,stimFlag] = gainModVRData(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% Returns mean deltaf averaged across trials and binned according to track
% length (Here: 2.5cm per bin)
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, numBin_max, nTracks, 2]
%       where numBin_max is the number of bins of the longest track
%       the fourth dimension is :
%           1 - results for gain trials
%           2 - results for non gain trials
%__________________________________________________________________________

stimFlag = [];


[ dataGainTrials , dataNonGainTrials, numBin_max ] = prepdataforfunction_RGroup('gainModVRData',signal, tt, auxlab, disp, vel, lagTime);
%  dataGainTrials and dataNonGainTrials are a cell arrays [numTracks, 1]
%       Each cell contains an array [nROIs, numBins, numTrials] where
%           numBins is the number of bins for the track
%           numTrials is the number of gain or other trials for the track
%       /!\ Tracks can be of different lengths, but here bins are fixed
%       (2.5cm) and hence numBins can be different between cells

% Number of tracks and ROIs
nTracks = length(dataGainTrials);
nROIs = size(signal,1);

% Initialise output array
X = nan(nROIs, numBin_max, nTracks, 2); 

for iTrack = 1:nTracks
    
    if ~isempty(dataGainTrials{iTrack})
        numBins = size(dataGainTrials{iTrack}, 2);
        % Average over trials
        X(:, 1:numBins, iTrack, 1) = nanmean(dataGainTrials{iTrack},3);
    end
    if ~isempty(dataNonGainTrials{iTrack})
        numBins = size(dataNonGainTrials{iTrack}, 2);
        % Average over trials
        X(:, 1:numBins, iTrack, 2) = nanmean(dataNonGainTrials{iTrack},3);
    end
end

end

function [X,stimFlag] = lmiVR(signal, tt, auxlab, ~, ~, ~, stimTimes, ustm)
% Computes the locomotion modulation index. Result is an unbiased estimate
% of the LMI over all stimuli.
VRdata = auxlab;

auxlab = VRdata.actlab;

exclude_stm = {'U'};
include_stm = {'none'};

X = compute_lmi(...
    signal, tt, auxlab, stimTimes, ustm, 0, [], include_stm, exclude_stm);

stimFlag = [];

end

function [X,stimFlag] = rwdOnsetData(signal, tt, auxlab, ~, ~, ~, ~, ~)
% Signal surrounding the onset of reward averaged over trials
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nSamples_RewardWindow, nRewardTypes, nTracks+1, 2]
%       For each ROI, returns the mean and sem of the average DF/F0 
%       surrounding the onset of reward (hit) or default reward, for each
%       track or average over all tracks (=>nTracks+1).
%       Dimensions: 
%           1 - ROI
%           2 - samples included in time window around reward
%           3 - reward type (reward or default reward)
%           4 - track number or all tracks (last)
%           5 - mean or sem
%
% Similar functions: RZdistOnsetData
%__________________________________________________________________________

stimFlag = [];

% 1) Convert function input to relevant Data 

[ dataDf_RewardWindow , window_samples, window_time, nTracks ] = prepdataforfunction_RGroup('rwdOnsetData',signal, tt, auxlab);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%   dataDf_RewardWindow : cell array [1, nRewardTypes]
%           Here for example nRewardTypes = 2
%               case 1: reward   case 2: default reward
%       A cell can be empty if there is no data for this type of reward.
%       Otherwise, each cell contains a cell array [1, nTracks], where
%           for each track the DF/F0 of each ROI was extracted for a time 
%           window (2*window_samples) anytime a reward event was detected, 
%           during all trials of the track (nSamples_event)
%           => array [nROIs, 2*window_samples, nSamples_event] 
%   window_samples*2 : number of samples included to extract the signal 
%       surrounding the onset of reward  
%   window_time : time corresponding to window_samples
%   nTracks : number of tracks
%       
% More Parameters
% ===============
f_resample = 40; % re-sampling frequency - Only works for 2 sec windows (40*4)
nSamples_RewardWindow = f_resample*4;

nRewardTypes = length(dataDf_RewardWindow);
nROIs = size(signal,1);

% Initialise output
X = nan(nROIs, nSamples_RewardWindow, nRewardTypes, nTracks+1, 2); 

for iRwdType = 1:nRewardTypes
    
    if isempty(dataDf_RewardWindow{iRwdType})
        continue
    end
    
    X_raw = nan(nROIs, window_samples*2, nTracks+1, 2);
    
    % Compute mean and sem for each track for each ROI 
    % (signal from onset of reward averaged over trials)
    for iTrack = 1:nTracks
        dataDf_track = dataDf_RewardWindow{iRwdType}{iTrack}; % size: [nROIs, 2*window_samples, nSamples_event] 
        if ~isempty(dataDf_track)
            X_raw(:,:,iTrack,1) = nanmean(dataDf_track,3);
            X_raw(:,:,iTrack,2) = nansem(dataDf_track,3);
        end
    end
    % Same, but for all tracks
    alltracks_dataDf = cat(3, dataDf_RewardWindow{iRwdType}{:}); % size: [nROIs, 2*window_samples, nSamples_event_alltracks] 
    X_raw(:,:,nTracks+1, 1) = nanmean(alltracks_dataDf ,3);
    X_raw(:,:,nTracks+1, 2) = nansem(alltracks_dataDf, 3);
    
    % Resample to 40Hz
    for iTrack = 1:nTracks+1
        for iData = 1:2
            
            if any(~isnan(X_raw(:,:,iTrack,iData)))
                
                data_resampled = ( resample(X_raw(:,:,iTrack,iData)',window_time,f_resample) )';
                nSamples_resamp = size(data_resampled,2);
                
                % Start filling resampled data 
                X(:, 1:min(nSamples_RewardWindow, nSamples_resamp), iRwdType, iTrack, iData) = data_resampled(:, 1:min(nSamples_RewardWindow, nSamples_resamp));
                
                % Check if data is missing
                if nSamples_resamp < nSamples_RewardWindow % pad with last datapoint
                    numPad = nSamples_RewardWindow - nSamples_resamp;
                    X(:, nSamples_resamp+1:end, iRwdType, iTrack, iData) = repmat(data_resampled(:,end), [1, numPad]);
                end 
            end 
        end  
    end
    
    % Replace first/last 5 data points as these can be weird from downsampling
    X(:,end-5:end,iRwdType,:,:) = X_raw(:,end-5:end,:,:);
    X(:,1:5,iRwdType,:,:) = X_raw(:,1:5,:,:);

end


end

function [X,stimFlag] = RZdistOnsetData(signal, tt, auxlab, disp, ~, ~,  ~, ~)
% Signal surrounding the onset of REWARD ZONE averaged over trials
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nSamples_RewardWindow, nRewardTypes, nTracks+1, 2]
%       For each ROI, returns the mean and sem of the average DF/F0 
%       surrounding the onset of reward (hit) or default reward, for each
%       track or average over all tracks (=>nTracks+1).
%       Dimensions: 
%           1 - ROI
%           2 - samples included in time window around reward
%           3 - reward type (reward or default reward)
%           4 - track number or all tracks (last)
%           5 - mean or sem
%
% Similar functions: rwdOnsetData
%__________________________________________________________________________

stimFlag = [];

% 1) Convert function input to relevant Data 

[ dataDf_RZoneWindow , window_samples, window_time, nTracks ] = prepdataforfunction_RGroup('rwdOnsetData',signal, tt, auxlab, disp);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%   dataDf_RZoneWindow : cell array [1, nRewardTypes]
%           Here for example nRewardTypes = 2
%               case 1: reward   case 2: default reward
%       A cell can be empty if there is no data for this type of reward.
%       Otherwise, each cell contains a cell array [1, nTracks], where
%           for each track the DF/F0 of each ROI was extracted for a time 
%           window (2*window_samples) around the reward zone anytime a 
%           reward event (type 1 or 2) was detected, considering all trials
%           of the track (nSamples_event)
%           => array [nROIs, 2*window_samples, nSamples_event] 
%   window_samples*2 : number of samples included to extract the signal 
%       surrounding the onset of reward  
%   window_time : time corresponding to window_samples
%   nTracks : number of tracks
%       
% More Parameters
% ===============
f_resample = 40; % re-sampling frequency - Only works for 2 sec windows (40*4)
nSamples_RewardWindow = f_resample*4;

nRewardTypes = length(dataDf_RZoneWindow);
nROIs = size(signal,1);

% Initialise output
X = nan(nROIs, nSamples_RewardWindow, nRewardTypes, nTracks+1, 2); 

for iRwdType = 1:nRewardTypes
    
    if isempty(dataDf_RZoneWindow{iRwdType})
        continue
    end
    
    X_raw = nan(nROIs, window_samples*2, nTracks+1, 2);
    
    % Compute mean and sem for each track for each ROI 
    % (signal from onset of reward averaged over trials)
    for iTrack = 1:nTracks
        dataDf_track = dataDf_RZoneWindow{iRwdType}{iTrack}; % size: [nROIs, 2*window_samples, nSamples_event] 
        if ~isempty(dataDf_track)
            X_raw(:,:,iTrack,1) = nanmean(dataDf_track,3);
            X_raw(:,:,iTrack,2) = nansem(dataDf_track,3);
        end
    end
    % Same, but for all tracks
    alltracks_dataDf = cat(3, dataDf_RZoneWindow{iRwdType}{:}); % size: [nROIs, 2*window_samples, nSamples_event_alltracks] 
    X_raw(:,:,nTracks+1, 1) = nanmean(alltracks_dataDf ,3);
    X_raw(:,:,nTracks+1, 2) = nansem(alltracks_dataDf, 3);
    
    % Resample to 40Hz
    for iTrack = 1:nTracks+1
        for iData = 1:2
            
            if any(~isnan(X_raw(:,:,iTrack,iData)))
                
                data_resampled = ( resample(X_raw(:,:,iTrack,iData)',window_time,f_resample) )';
                nSamples_resamp = size(data_resampled,2);
                
                % Start filling resampled data 
                X(:, 1:min(nSamples_RewardWindow, nSamples_resamp), iRwdType, iTrack, iData) = data_resampled(:, 1:min(nSamples_RewardWindow, nSamples_resamp));
                
                % Check if data is missing
                if nSamples_resamp < nSamples_RewardWindow % pad with last datapoint
                    numPad = nSamples_RewardWindow - nSamples_resamp;
                    X(:, nSamples_resamp+1:end, iRwdType, iTrack, iData) = repmat(data_resampled(:,end), [1, numPad]);
                end 
            end 
        end  
    end
    
    % Replace first/last 5 data points as these can be weird from downsampling
    X(:,end-5:end,iRwdType,:,:) = X_raw(:,end-5:end,:,:);
    X(:,1:5,iRwdType,:,:) = X_raw(:,1:5,:,:);

end


end

function [X,stimFlag] = successRateSMI_VR(signal, ~, auxlab, disp, ~, ~, ~, ~)
% Description 
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nTracks]
%       SMI: spatial modulation index. Calculated as percent success rate over trials
%       divided by the shuffled sucess rate for each VR track
%       Note that all rows are just copies as success rate applies to all ROIs.
%
% see also: Shuffle
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 

[ success_rate, dataReward, dataLicks, dataDistance, TracksParameters ] = prepdataforfunction_RGroup('successRateSMI_VR', auxlab, disp);

% 2) Analysis 
%
% Format of Data for analysis
% ===========================
%   success_rate : array [1, numTracks]
%           percent success rate over trials for the track
%           NOTE: if there is no data for the track, the entry is NaN
%   dataReward : cell array [1, numTracks]
%       Each cell contains a cell array [1, numTrials] where numTrials is 
%       the number of trials of the track
%   dataLicks : cell array [1, numTracks]
%       Same as dataReward
%   dataDistance : cell array [1, numTracks]
%       Same as dataReward
%   TracksParameters is a structure array [1, numTracks] with fields:
%       - rewardZone (integer)
%       - rewardDefault (integer)
%       These two are parameters to identify the reward zone where 
%       the animal is supposed to lick (see below)
%
% NOTE ABOUT REWARD DATA: 
%           0 = no reward   
%           1 = reward  
%           2 = default reward
%       
% Parameters for analysis
% =======================
numShuffles = 1000;

numTracks = length(TracksParameters);
shuff_rate = nan(numShuffles, numTracks);

for iTrack = 1:numTracks

    if isnan(success_rate(iTrack))
        continue
    end
    
    % Get params of track
    rz = TracksParameters(iTrack).rewardZone;
    default = TracksParameters(iTrack).rewardDefault;
    numTrials = length(dataReward{iTrack});
    
    hit_shuffle = zeros(numShuffles, numTrials);
    
    for iTrial = 1:numTrials
        
        trial_length = length(dataReward{iTrack}{iTrial});
        
        if any(dataReward{iTrack}{iTrial}==1) % hit trial
            Idx_startRwd = find(dataReward{iTrack}{iTrial}==1,1);
            
        elseif any(dataReward{iTrack}{iTrial}==2) % default trial
            Idx_startRwd = find(dataReward{iTrack}{iTrial}==2,1);
            
        else % No reward
            Idx_startRwd = trial_length;
        end
        
        prerwdLicks = dataLicks{iTrack}{iTrial}(1:Idx_startRwd);
        pad = zeros(1,trial_length-Idx_startRwd);
        tmpdata = [pad, prerwdLicks];
        
        shufflick = Shuffle(repmat(tmpdata,numShuffles,1),2);
        
        selected_track_zone = dataDistance{iTrack}{iTrial}>=rz & dataDistance{iTrack}{iTrial}<default ;
        hit_shuffle(:,iTrial) = double( sum(shufflick(:, selected_track_zone), 2)>0 );
    end
    
    shuff_rate(:,iTrack) = sum(hit_shuffle,2)/numTrials; % percent success rate, num hits, num misses
end

% Average over numShuffles
shuff_rate = nanmean(shuff_rate,1); % percent success rate, num hits, num misses

% Ratio
R = success_rate./shuff_rate;

% Simple copy for Output of main function
nROIs = size(signal,1);
X = repmat(R, nROIs, 1);

end

function [X,stimFlag] = pairedTestGratVSrwdZone(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% Description
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nTracks, 3]
%       Third dimension: For each ROI for each track, 
%           1 - test pre- versus post- reward zone
%               (Wilcoxon signed rank test for zero median)
%           2 - average DF/F0 pre-reward zone
%           3 - average DF/F0 post-reward zone
%
% see also: signrank
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 
   
[RZdatapre, RZdatapost] = prepdataforfunction_RGroup('pairedTestGratVSrwdZone',signal, tt, auxlab, disp, vel, lagTime);

% 2) Analysis 
%
% Format of Data for analysis
% ===========================
%  Data should be two cell arrays [numTracks, 1]
%       - RZdatapre : DF/F0 for pre-reward zone
%       - RZdatapost : DF/F0 for post-reward zone
%       Each cell contains an array [nROIs, numTrials] where numTrials is
%       the number of trials for the track. The values DF/F0 were obtained
%       by averaging, for each trial, over selected bins pre- / post- rwd.
%       
% Parameters for analysis
% =======================
alpha = 0.001; % significance level for p-value

% Number of tracks and ROIs
nTracks = size(RZdatapre,1);
nROIs = size(signal,1);

% Initialise output array
X = nan(nROIs, nTracks, 3); 

for iTrack = 1:nTracks
    
    if ~isempty(RZdatapre{iTrack})
        
        for iROI = 1:nROIs
            X(iROI,iTrack,1) = signrank(RZdatapre{iTrack}(iROI,:),RZdatapost{iTrack}(iROI,:), 'alpha', alpha);
        end
        
        X(:,iTrack,2) = nanmean(RZdatapre{iTrack},2);
        X(:,iTrack,3) = nanmean(RZdatapost{iTrack},2);
        
    end
    

end

end

function [X,stimFlag] = tempMatchVR_HitMiss(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% Description
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 2, nTracks]
%       Note that all rows are just copies as decoder applies to all ROIs.
%       First column: For each track, give the percentage of correctly 
%           classified trials, either Hit or Miss, given the activities.
%       Second column: Chance level.
%
% see also: templatematching
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 

[dataDfHit, dataDfMiss] = prepdataforfunction_RGroup('tempMatchVR_HitMiss',signal, tt, auxlab, disp, vel, lagTime);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%  Data should be two cell arrays [numTracks, 1]
%       - dataDfHit : average DF/F0 for hit trials
%       - dataDfMiss : average DF/F0 for miss trials
%       Each cell contains an array [nROIs, numTrials] where numTrials is
%       the number of hit or miss trials.
     

% Number of tracks and ROIs
nTracks = size(dataDfHit,1);
nROIs = size(signal,1);

% Initialise output
X = nan(1, 2, nTracks); 

for iTrack = 1:nTracks
    
    if ~isempty(dataDfHit{iTrack})
        
        numTrials_max = max(size(dataDfHit{iTrack},2), size(dataDfMiss{iTrack},2));
        
        dataDfPad = nan(nROIs, 2, numTrials_max); 
                
        for iROI = 1:nROIs
            dataTmp = {dataDfHit{iTrack}(iROI,:) , dataDfMiss{iTrack}(iROI,:)};
            dataDfPad(iROI, :, :) = padcat(dataTmp{:});  
        end
        
        X(1,1,iTrack) = templatematching(dataDfPad); % decoder accuracy
        X(1,2,iTrack) = 1/2; % chance level of decoder (i.e. determined by how many track inputs)
        
    end
    
end

X = repmat(X, nROIs, 1, 1);

end

function [X,stimFlag] = tempMatchVR_RZ_CuedUncued(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% Description
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 2, nTrackTypes]
%       Note that all rows are just copies as decoder applies to all ROIs.
%       First column: For each track pair, give the percentage of correctly 
%           classified trials (e.g. visually cued vs uncued tracks).
%       Second column: Chance level.
%
% see also: templatematching
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 

[ RZdataDf, trackType_Idx ] = prepdataforfunction_RGroup('tempMatchVR_RZ_CuedUncued',signal, tt, auxlab, disp, vel, lagTime);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%   RZdataDf is a cell array [numTracks, 1] 
%       Each cell contains an array [nROIs, numTrials] with average DF/F0 
%       for hit trials, where numTrials is the number of hit trials.
%   trackType_Idx is a cell array [nTrackTypes, 1]
%       For each track type, gives the index of corresponding tracks

% Number of tracks and ROIs
nTracks = length(RZdataDf);
nROIs = size(signal,1);
nTrials = cellfun(@(x) size(x,2), RZdataDf);
nTrials_max = max(nTrials);
nTrackTypes = length(trackType_Idx);

RZdataDfPad = nan(nROIs, nTracks, nTrials_max);
for iTrack = 1:nTracks
    if ~isempty(RZdataDf{iTrack})
        RZdataDfPad(:, iTrack, 1:nTrials(iTrack)) = RZdataDf{iTrack};
    end
end
        
% Initialise output
X = nan(1, 2, nTrackTypes); 

for iTrackType = 1:nTrackTypes
    
    RZdataDfPad_currentType = RZdataDfPad(:,trackType_Idx{iTrackType},:);
    % discard tracks that are all nan
    keep_track = any(any(~isnan(RZdataDfPad_currentType),3),1);
    RZdataDfPad_currentType = RZdataDfPad_currentType(:,keep_track,:);
    nTracks_currentType = size(RZdataDfPad_currentType,2);
    if nTracks_currentType>1
        X(1,1,iTrackType) = templatematching(RZdataDfPad_currentType); % decoder accuracy
        X(1,2,iTrackType) = 1/nTracks_currentType; % chance level of decoder (i.e. determined by how many track inputs)
    end
    
end

X = repmat(X, nROIs, 1, 1);

end

function [X,stimFlag] = trialTime(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% Returns trial time statistics
% (median/mean/sem/min/max/numTrials) by tracks matrix
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 6, nTracks]
%       Note that all rows are just copies as the stats apply to all ROIs.
%__________________________________________________________________________


% input handeling
VRlab = auxlab;
% Get trial times
[~, ~, ~, ~, ~, ~, ~, ~, trialTimes] = VRdata2tracks(signal, tt, VRlab, disp, vel, lagTime);

nTracks =  length(trialTimes);
nROIs = size(signal,1);

X = nan(1, 6, nTracks);
for iTrack = 1:nTracks
    if ~isempty(trialTimes{iTrack})
        X(1, 1, iTrack) = nanmedian(trialTimes{iTrack});
        X(1, 2, iTrack) = nanmean(trialTimes{iTrack});
        X(1, 3, iTrack) = nansem(trialTimes{iTrack});
        X(1, 4, iTrack) = min(trialTimes{iTrack});
        X(1, 5, iTrack) = max(trialTimes{iTrack});
        X(1, 6, iTrack) = size(trialTimes{iTrack},1);
    end
end
X = repmat(X, nROIs, 1, 1); % data is same for each ROI

stimFlag = [];

end

function [X,stimFlag] = slowMedfastTrials(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% Description
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nTracks, 6]
%       third dimension
%           1 - pval Wilcoxon rank sum test
%           2 - mean DF/F0 for slow trials 
%           3 - mean DF/F0 for med trials 
%           4 - mean DF/F0 for fast trials 
%           5 - trial time for 75th percentile (Slow trials) 
%           6 - trial time for 25th percentile (Fast trials)
%__________________________________________________________________________


stimFlag = []; 

% 1) Convert function input to relevant Data 

[trialTimes, dataDf_RZone_slow, dataDf_RZone_med, dataDf_RZone_fast] = prepdataforfunction_RGroup('slowMedfastTrials',signal, tt, auxlab, disp, vel, lagTime);


% 2) Analysis
%
% Format of Data for analysis
% ===========================
%   trialTimes is a cell array [nTracks, 1] 
%       Each cell contains an array [numTrials, 1]
%   dataDf_RZone_X is a cell array [nTracks, 1] 
%       Each cell contains an array [nROIs, numTrials] with average DF/F0 
%       surrounding reward zone for slow, med, or fast trials.

nTracks = size(dFTrials,1);
nROIs = size(signal,1);

X = nan(nROIs, nTracks, 6); 

if sum(~cellfun(@isempty, dataDf_RZone_fast))>0
    dataTmp_slow = cell(1,nTracks); 
    dataTmp_med = cell(1,nTracks);
    dataTmp_fast = cell(1,nTracks);

    for iROI = 1:nROIs
        for iTrack = 1:nTracks
            dataTmp_slow{iTrack} = dataDf_RZone_slow{iTrack}(iROI,:);
            dataTmp_med{iTrack} = dataDf_RZone_med{iTrack}(iROI,:);
            dataTmp_fast{iTrack} = dataDf_RZone_fast{iTrack}(iROI,:);
        end
        RZdataPad_slow(:,:,iROI) = padcat(dataTmp_slow{:});  %#ok<AGROW>
        RZdataPad_med(:,:,iROI) = padcat(dataTmp_med{:});  %#ok<AGROW>
        RZdataPad_fast(:,:,iROI) = padcat(dataTmp_fast{:});  %#ok<AGROW>
    end
    RZdataPad_slow = permute(RZdataPad_slow, [3,1,2]);
    RZdataPad_med = permute(RZdataPad_med, [3,1,2]);
    RZdataPad_fast = permute(RZdataPad_fast, [3,1,2]);
    
    for iTrack = 1:nTracks
        if sum(~isnan(squeeze(RZdataPad_slow(iROI,iTrack,:))))>0
            
            for iROI = 1:nROIs
                dataTable = [...
                    squeeze(RZdataPad_slow(iROI,iTrack,:)) ; ...
                    squeeze(RZdataPad_fast(iROI,iTrack,:)) ; ...
                    squeeze(RZdataPad_med(iROI,iTrack,:)) ];
                
                group = [...
                    ones( size(RZdataPad_slow(iROI,iTrack,:),3) ,1) ;...
                    2*ones( size(RZdataPad_fast(iROI,iTrack,:),3) ,1) ; ...
                    3*ones( size(RZdataPad_med(iROI,iTrack,:),3) ,1) ];
                
                X(iROI,iTrack,1) = anova1(dataTable, group, 'off'); % Wilcoxon rank sum test
                X(iROI,iTrack,2) = nanmean(squeeze(RZdataPad_slow(iROI,iTrack,:))); % Slow Trials
                X(iROI,iTrack,3) = nanmean(squeeze(RZdataPad_med(iROI,iTrack,:))); % Med trials
                X(iROI,iTrack,4) = nanmean(squeeze(RZdataPad_fast(iROI,iTrack,:))); % Fast trials
                X(iROI,iTrack,5) = prctile(trialTimes{iTrack}, 75); % trial time for 75th percentile (Slow trials)
                X(iROI,iTrack,6) = prctile(trialTimes{iTrack}, 25); % trial time for 25th percentile (Fast trials)
            end
        end
    end
    
end

end

function [X,stimFlag] = mancovaVR(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% see: mancovan 
% [ T, p, FANCOVAN, pANCOVAN, stats ] = mancovan(Y, groups, covariates, options)
%     Y          - [ N x M ] (double) - multivariate response
%     groups     - [ N x G ] (int)    - qualitative variables
%     covariates - [ N x C ] (double) - quantitative variables
% Options:
%     'group-group'         - include group-group interactions
%     'covariate-covariate' - include covariate-covariate interactions
%     'group-covariate'     - include group-covariate interactions
%     'over-determined'     - use over-determined coding for the design matrix
%     'sigma-restricted'    - use sigma-restricted coding for the design matrix
%     'SVD'                 - reduce the dimensionality of Y using an SVD
%     'verbose'             - display extra information to the command window
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 

[trialTimes, hitTrial, dataDf_RZ, dataVel_RZ] = prepdataforfunction_RGroup('mancova',signal, tt, auxlab, disp, vel, lagTime);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%   trialTimes is a cell array [nTracks, 1] 
%       Each cell contains an array [numTrials, 1]
%   hitTrial is a cell array [nTracks, 1] 
%       Declares for each trial of each track if the animal got a reward.
%       Each cell contains a boolean vector [numTrials, 1]
%   dataDf_RZ is a cell array [nTracks, 1] 
%       Each cell contains an array [nROIs, numTrials] with average DF/F0 
%       surrounding reward zone for slow, med, or fast trials.


% Get number of tracks and ROIs
nTracks = size(dFTrials,1);
nROIs = size(signal,1);

X = nan(nROIs, 4, nTracks); 
Y = nan(nROIs, 4, nTracks); 
Z = nan(nROIs, 4, nTracks); 

for iTrack = 1:nTracks
    num_hit_trials = sum(hitTrial{iTrack});
    if ~isempty(dataDf_RZ{iTrack}) && num_hit_trials~=size(dataDf_RZ{iTrack},2) && num_hit_trials>=3
        
        data = dataDf_RZ{iTrack}'; % [numTrials, nROIs]
        
        groups = double(hitTrial{iTrack}); % group 1 hit - size [trials x 1]
        groups(groups==0)=2; % group 2 miss - size [trials x 1]
        
        covariates = [dataVel_RZ{iTrack}, trialTimes{iTrack}]; 
        
        [ ~, Z(1,1,iTrack), ~, ~, ~] = mancovan(data, groups, [], {'SVD'});

        [ ~, Z(1,2:4,iTrack)] = mancovan(data, groups, covariates, {'SVD'});

        for iROI = 1:nROIs
            data = dataDf_RZ{iTrack}(iROI,:)'; % [numTrials, nROIs]
            
            [ Y(iROI,1,iTrack), X(iROI,1,iTrack)] = mancovan(data, groups);
            
            [ TgpCo, pgpCo] = mancovan(data, groups, covariates);
            Y(iROI,2:4,iTrack) = TgpCo';
            X(iROI,2:4,iTrack) = pgpCo';

        end
    end
    Z(:,:,iTrack) = repmat(Z(1,:,iTrack), nROIs, 1, 1);
end

X = [X, Z];


end

function [X,stimFlag] = mancovaVRshuffDist(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% see: mancovan 
% [ T, p, FANCOVAN, pANCOVAN, stats ] = mancovan(Y, groups, covariates, options)
%     Y          - [ N x M ] (double) - multivariate response
%     groups     - [ N x G ] (int)    - qualitative variables
%     covariates - [ N x C ] (double) - quantitative variables
% Options:
%     'group-group'         - include group-group interactions
%     'covariate-covariate' - include covariate-covariate interactions
%     'group-covariate'     - include group-covariate interactions
%     'over-determined'     - use over-determined coding for the design matrix
%     'sigma-restricted'    - use sigma-restricted coding for the design matrix
%     'SVD'                 - reduce the dimensionality of Y using an SVD
%     'verbose'             - display extra information to the command window
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 

is_shuffle = true;

[trialTimes, hitTrial, dFTrials_perm, dataVel_RZ] = prepdataforfunction_RGroup('mancova',signal, tt, auxlab, disp, vel, lagTime, is_shuffle);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%   trialTimes is a cell array [nTracks, 1] 
%       Each cell contains an array [numTrials, 1]
%   hitTrial is a cell array [nTracks, 1] 
%       Declares for each trial of each track if the animal got a reward.
%       Each cell contains a boolean vector [numTrials, 1]
%   dFTrials_perm is a cell array [nTracks, 1] 
%       Each cell contains an array [nROIs, numTrials, numBinSamples] with DF/F0 

num_iters = 100;

% Get number of tracks and ROIs
nTracks = length(dFTrials_perm);
nROIs = size(signal,1);

X = nan(nROIs, 4, nTracks, num_iters); 
Z = nan(nROIs, 4, nTracks, num_iters); 

for iIter = 1:num_iters
    % Get random 8 bins from virtual corridor (first 30 bins) instead of
    % reward zone
    dataRange = randi(30, [1,8]);

    dataDf_RZ = cellfun(@(x) nanmean(x(:,:,dataRange),3), dFTrials_perm, 'UniformOutput', false);

    for iTrack = 1:nTracks
        
        num_hit_trials = sum(hitTrial{iTrack});
        
        if ~isempty(dataDf_RZ{iTrack}) && num_hit_trials~=size(dataDf_RZ{iTrack},2) && num_hit_trials>=3
            
            data = dataDf_RZ{iTrack}'; % [numTrials, nROIs]
            data(isnan(data)) = 0;
            
            groups = double(hitTrial{iTrack}); % group 1 hit - size [trials x 1]
            groups(groups==0)=2; % group 2 miss
            
            covariates = [dataVel_RZ{iTrack}, trialTimes{iTrack}]; 
            covariates(isnan(covariates)) = 0;
            [ ~, Z(1,1,iTrack,iIter)] = mancovan(data, groups, [], {'SVD'}); % , FANCOVAN, pANCOVAN, stats

            [ ~, Z(1,2:4,iTrack,iIter)] = mancovan(data, groups, covariates, {'SVD'});
         
            for iROI = 1:nROIs
                data = dataDf_RZ{iTrack}(iROI,:)'; % [numTrials, nROIs]
                data(isnan(data)) = 0;
                [ ~, X(iROI,1,iTrack,iIter)] = mancovan(data, groups);
         
                try
                    [ ~, pgpCo] = mancovan(data, groups, covariates);
                catch ME
                    if strcmp(ME.identifier, 'MATLAB:betainc:XOutOfRange')
                        dataRange = randi(30, [1,8]); % from end of trial; RZ starts 16 bins from end;
                        dataDf_RZ = cellfun(@(x) nanmean(x(:,:,dataRange),3), dFTrials_perm, 'UniformOutput', false);
                        
                        data = dataDf_RZ{iTrack}(iROI,:)';
                        data(isnan(data)) = 0;
                        [ ~, X(iROI,1,iTrack,iIter)] = mancovan(data, groups);
                        [ ~, pgpCo] = mancovan(data, groups, covariates);
                    end
                end
                X(iROI,2:4,iTrack,iIter) = pgpCo';
            
            end
            
        end
        Z(:,:,iTrack,:) = repmat(Z(1,:,iTrack,:), nROIs,1,1,1);
    end
    
end

X = [X, Z];

end

function [X,stimFlag] = peakVarByOnset(signal, tt, auxlab, disp, ~, ~, ~, ~)

VRdata = auxlab;
smp_dt_align = mean(diff(tt));
onsetTypes = {'distance', 'reward', 'start'} ;

rwdZoneSize = 40;
tWINDOW = 30;
trackNum = [1,2,3,5,7,8];

nRois = size(signal,1);
nOnsetTypes = length(onsetTypes);
nTracks = length(trackNum);

% Init output
X = nan(nRois, nOnsetTypes, nTracks);
stimFlag = [];

numTrials = max(VRdata.trialNum);
maxTrialLgth = 0;
for iTrial = 1:numTrials
    is_current_trial = VRdata.trialNum == iTrial;
    maxTrialLgth = max(maxTrialLgth, sum(is_current_trial));
end

is_reward = VRdata.reward ~= 0; % success and miss trials together
% NOTE ABOUT REWARD DATA:
%           0 = no reward
%           1 = reward (success trial)
%           2 = default reward (miss trial)

for iOnset = 1:length(onsetTypes)
    
    onsetType = onsetTypes{iOnset};
    
    % First get the relevant time window
    trialTime = nan(numTrials, maxTrialLgth);
    onsetCol = nan(numTrials,1);
    startCol = nan(numTrials,1);
    endCol = nan(numTrials,1);
    is_dataTrialType_reward = false(numTrials,1);
    
    for iTrial = 1:numTrials
        
        is_current_trial = VRdata.trialNum == iTrial;
        if ~any(is_current_trial)
            continue
        end
        
        if any( is_reward(is_current_trial) )
            is_dataTrialType_reward(iTrial) = true;
        end
        
        % Normalize desired onset time to 0
        if strcmp(onsetType, 'start')
            trialOnset = find(is_current_trial,1);
            trialTime(iTrial, 1:sum(is_current_trial)) = tt(1,is_current_trial) - tt(1,trialOnset);
            
            onsetCol(iTrial) = knnsearch(trialTime(iTrial,:)',0);
            startCol(iTrial) = knnsearch(trialTime(iTrial,:)',0);
            endCol(iTrial) = knnsearch(trialTime(iTrial,:)',tWINDOW*2);
            
        else
            if strcmp(onsetType, 'reward')
                rwdOnset = mean(tt(is_current_trial & is_reward));
                trialTime(iTrial, 1:sum(is_current_trial)) = tt(1,is_current_trial) - rwdOnset;
                
            else % strcmp(onsetType, 'distance')
                rwdZoneStart = (round(max(disp(is_current_trial))/10)*10) - rwdZoneSize;
                displmntTmp = disp;
                displmntTmp(~is_current_trial) = nan;
                distOnset = knnsearch(displmntTmp', rwdZoneStart);
                trialTime(iTrial, 1:sum(is_current_trial)) = tt(1,is_current_trial) - tt(1,distOnset);                
            end
            
            onsetCol(iTrial) = knnsearch(trialTime(iTrial,:)',0);
            startCol(iTrial) = knnsearch(trialTime(iTrial,:)',-tWINDOW);
            endCol(iTrial) = knnsearch(trialTime(iTrial,:)',tWINDOW);
            
            onsetCol(onsetCol ==1) = nan;
        end
    end
    
    maxAlignTrialLgth = max(diff([startCol, endCol]'));
    if mod(maxAlignTrialLgth,2) ~= 0
        maxAlignTrialLgth = maxAlignTrialLgth+1;
    end
    if strcmp(onsetType, 'start')
        tstart = onsetCol;
        tend = onsetCol + maxAlignTrialLgth;
    else
        tstart = onsetCol - maxAlignTrialLgth/2;
        tend = onsetCol + maxAlignTrialLgth/2;
    end
    
    
    
    % Now loop over track types to fill X
    for itrackNum = 1:nTracks
        
        is_current_track = VRdata.trackNum == trackNum(itrackNum);
        if ~any(is_current_track)
            continue
        end
        
        if strcmp(onsetType, 'start')
            ttTemplate = 0:smp_dt_align:tWINDOW*2;
        else
            ttTemplate = -tWINDOW:smp_dt_align:tWINDOW;
        end
        while length(ttTemplate) < maxAlignTrialLgth
            ttTemplate = [ttTemplate, ttTemplate(end)+smp_dt_align]; %#ok<AGROW>
        end
        while length(ttTemplate) > maxAlignTrialLgth
            ttTemplate = ttTemplate(1:end-1);
        end
        
        for iROI = 1:nRois
            % Separate trials
            alignSignal = nan(numTrials, maxAlignTrialLgth);
            alignTime = nan(numTrials, maxAlignTrialLgth);
            
            for iTrial = 1:numTrials
                if isnan(tstart(iTrial))
                    continue
                end
                is_current_trial = VRdata.trialNum == iTrial;
                is_current_trial_of_track = is_current_trial & is_current_track;
                if ~any(is_current_trial_of_track)
                    continue
                end
                
                trialData = signal(iROI, is_current_trial_of_track);
                
                if tend(iTrial)>maxTrialLgth
                    alignSignal(iTrial,:) = [trialData(tstart(iTrial)+1:maxTrialLgth), nan(1,abs(tend(iTrial)-maxTrialLgth))];
                    alignTime(iTrial,:) = [trialTime(iTrial, tstart(iTrial)+1:maxTrialLgth), nan(1,abs(tend(iTrial)-maxTrialLgth))];
                    
                elseif tstart(iTrial)>=1
                    alignSignal(iTrial,:) = trialData(tstart(iTrial):tend(iTrial)-1);
                    alignTime(iTrial,:) = trialTime(iTrial, tstart(iTrial):tend(iTrial)-1);
                    
                elseif tstart(iTrial)==0
                    alignSignal(iTrial,:) = trialData(tstart(iTrial)+1:tend(iTrial));
                    alignTime(iTrial,:) = trialTime(iTrial, tstart(iTrial)+1:tend(iTrial));
                    
                else
                    alignSignal(iTrial,:) = [nan(1,abs(tstart(iTrial))+1), trialData(1:tend(iTrial)-1)];
                    alignTime(iTrial,:) = [nan(1,abs(tstart(iTrial))+1), trialTime(iTrial, 1:tend(iTrial)-1)];
                end
            end
            dataTable = alignSignal(is_dataTrialType_reward,:);
            dataTable(sum(isnan(dataTable), 2)==size(dataTable,2),:)=[];
            
            [~,maxPeakIdx] = max(dataTable,[],2);
            peakTimes = ttTemplate(maxPeakIdx);
            
            X(iROI,iOnset,itrackNum) = var(peakTimes);
            
        end
    end
end

end