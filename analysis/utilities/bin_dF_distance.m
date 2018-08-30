% Compute signals for each ROI for all trials of one track.
% 
% Calculate mean dF/F0 for each bin (distance or time bin)
% To map dF/F samples to location or track time, each trial is divided into 
% numBins bins, and the mean of all samples within that bin is calculated.
% The standard error of the mean (SEM) represents the SEM of mean dF 
% between trials, NOT within trials
%
% Notations
% ---------
%   nSamples = Total number of samples for the track (all trials)
%   numTotTrials = Total number of trials for the track 
%        
% Inputs 
% ------
%   VRdata : array [9, nSamples] with rows  / or  structure fields:
%       Each row/field is a vector of size [1, nSamples]
%       1 - time : Time of samples since the start of the recording (s)
%       2 - trackPosition : Distance travelled since the start of the 
%                           recording (cm)
%       3 - frameLagTime : Instantaneous frame by frame lag in VR (s)
%       4 - velocity : Instantaneous velocity (cm/s)
%       5 - trackNum
%       6 - reward : 0 = no reward ; 1 = got reward ; 2 = missed reward
%       7 - trialNum
%       8 - licks : 0 = no lick ; 1+ = licks 
%       9 - actlab
%   dFdata : array [nROIs, nSamples]
%       Fluorescence activity traces per ROI for each trial for the track
%   numBins :  integer   (e.g. 48, 64, or 96)
%       Number of bins for the current track
%       Tracks can be of length 120, 160 or 240cm; Bins are 2.5cm
%   binBy : string   (Default: 'distance')   
%       Possible values = 'distance' , 'time'
%        
% Inputs (optional)
% ----------------
%   dataTrialsOnly : boolean   (Default: false)
%       Return only trials with data (remove nan trials of other tracks or
%       categories)
%   useHitsOnly : boolean   (Default: false)
%       Some descript
%   gainTrials : boolean   (Default: false)
%       Some descript
%   nonGainTrials : boolean   (Default: false)
%       Some descript
%   useMissOnly : boolean   (Default: false)
%       Some descript
%
% Outputs
% -------
%   dFTrials : array [nROIs, numBins, numTotTrials]
%       Binned fluorescence dF/F signal for each trial of the track
%   meandF : array [nROIs, numBins]
%       Average dF/F over trials
%   semdF : array [nROIs, numBins]
%       Standard error of the mean 'meandF' over trials
%   hitTrial : Boolean vector [numTotTrials, 1]
%       Declares for each trial of each track if the animal got a reward.
%   velocityTrials : array [1, numBins, numTotTrials]
%       Binned velocity data 
%   trialTimes : array [1, numTotTrials]
%       Get total length of each trial, not binned by track length 
%   rwdTrial : array [1, numBins, numTotTrials]
%       Binned reward data 
%   licksTrial : array [1, numBins, numTotTrials]
%       Binned licks data 
%   actlabTrial : array [1, numBins, numTotTrials]
%       Binned action label data (related to velocity)
%
% See also: VRdata2tracks


function [dFTrials, meandF, semdF, hitTrial, velocityTrials,...
    trialTimes, rwdTrial, licksTrial, actlabTrial] = ...
    bin_dF_distance(VRdata, dFdata, numBins, binBy, ...
        dataTrialsOnly, useHitsOnly, gainTrials, nonGainTrials, useMissOnly)

% ** Check optional input
if nargin<9 || isempty(useMissOnly)
    useMissOnly = false;
end
if nargin<8 || isempty(nonGainTrials)
    nonGainTrials = false;
end
if nargin<7 || isempty(gainTrials)
    gainTrials = false;
end
if nargin<6 || isempty(useHitsOnly)
    useHitsOnly = false;
end
if nargin<5 || isempty(dataTrialsOnly)
    dataTrialsOnly = false;
end
% ** End check optional input

if nonGainTrials || gainTrials
    nthTrial = 5; % e.g. gain trials are every 5th trial starting with trial 4 (trialStart value below).
    trialStart = 4; 
end

% If VRdata input is a structure, convert into array 
if isstruct(VRdata)
    VRdata = struct2cell(VRdata);
    VRdata = cell2mat(VRdata);
end

% Decide how to bin data
if strcmp(binBy, 'time')
    binData = VRdata(1,:);
elseif strcmp(binBy, 'distance')
    binData = VRdata(2,:);
else
    error('Criteria to bin by is not recognized')
end

trialsNumber = VRdata(7,:);
list_trialNumbers = unique(trialsNumber);
numTotTrials = length(list_trialNumbers);
rewardData =VRdata(6,:); 
velocityData = VRdata(4,:);
licksData = VRdata(8,:);
actlabData = VRdata(9,:);
ttData = VRdata(1,:);
nROIs = size(dFdata,1);

% Initialize output data - organised by trials
dFTrials = nan(nROIs, numBins, numTotTrials);
hitTrial = false(numTotTrials,1);
dataTrial = ones(numTotTrials,1); % 1 = trial has data ; 0 = does not
velocityTrials = nan(1, numBins, numTotTrials);
rwdTrial = nan(size(rewardData,1), binnr, numTotTrials);
licksTrial = nan(size(licksData,1), binnr, numTotTrials);
actlabTrial = nan(size(actlabData,1), binnr, numTotTrials);
trialTimes = nan(numTotTrials,1);

% Loop over trials of the track
for iTrial = 1:numTotTrials
    
    selected_trial = list_trialNumbers(iTrial);
    
    % Collect data of trial
    trialIdx = trialsNumber == selected_trial; % vector with boolean elements
    if ~any(trialIdx) % 'any' finds if there is at least 1 true element
        dataTrial(iTrial,1) = 0; 
        continue;
    end
    
     % get trial times (total length of each trial, not binned by track length)
    tttmp = ttData(:,trialIdx);
    trialTimes(iTrial) = tttmp(end)-tttmp(1);
    
    % Get all other data binned over track length
    trialSignal = dFdata(:,trialIdx); % [nROIs, numDataPts_trial_k]
    numDataPts_trial_k = size(trialSignal,2);
    trialspeed = velocityData(1,trialIdx);
    trialrwd = rewardData(:,trialIdx);
    triallick = licksData(:,trialIdx);
    trialactlab = actlabData(:,trialIdx);
    binInfo = binData(1,trialIdx);
    
    [~,bin] = histc(binInfo,linspace(min(binInfo),max(binInfo),numBins));
    % Size of bin is [1, numBins]
    
    % Bin deltaf data for each ROI
    r_i = 1:numDataPts_trial_k;
    for iROI = 1:nROIs
        ty = sparse(r_i, bin, trialSignal(iROI,:)); % size: [numDataPts_trial_k, numBins]
        dFTrials(iROI,:,iTrial) = full(sum(ty,1)./sum(ty~=0,1)); % size: [1, numBins]
    end
    
    % Bin velocity data
    ty = sparse(r_i, bin, trialspeed(1,:));
    velocityTrials(1,:,iTrial) = full(sum(ty,1)./sum(ty~=0,1));
    % bin rwdard data
    ty = sparse(1:length(binInfo),bin,trialrwd(1,:));
    rwdTrial(1,:,iTrial) = full(sum(ty)./sum(ty~=0));
    % bin licks data
    ty = sparse(1:length(binInfo),bin,triallick(1,:));
    licksTrial(1,:,iTrial) = full(sum(ty)./sum(ty~=0));
    % bin actlab data
    ty = sparse(1:length(binInfo),bin,trialactlab(1,:));
    actlabTrial(1,:,iTrial) = full(sum(ty)./sum(ty~=0));
    
    % Check if animal got reward during trial
    hitTrial(iTrial) = any(rewardData(trialIdx)==1);
    
end

% if nan values (i.e. no licking on a certain trial), make ==0 for reward and licks
rwdTrial(isnan(rwdTrial))=0;
licksTrial(isnan(licksTrial))=0;

% Return data based on Hit vs Miss or gain vs nongain outcome
if useHitsOnly && useMissOnly
    error('Can not request data from only hits and misses at the same time... pick one (or neither [both set to zero] to get all data)')
end
if gainTrials && nonGainTrials
    error('Can not request data from only gain trial and non-gain trials at the same time... pick one (or neither [both set to zero] to get all data)')
end

% first pick out gain trials (or non-gain trials)
if gainTrials
    gainTrialsIdx = false(1,size(dFTrials,3));
    gainTrialsIdx(1,trialStart:nthTrial:size(dFTrials,3)) = true;
    dFTrials = dFTrials(:,:,gainTrialsIdx);
    velocityTrials = velocityTrials(:,:,gainTrialsIdx);
    licksTrial = licksTrial(:,:,gainTrialsIdx);
    rwdTrial = rwdTrial(:,:,gainTrialsIdx);
    actlabTrial = actlabTrial(:,:,gainTrialsIdx);
    trialTimes = trialTimes(gainTrialsIdx,1);
    dataTrial = dataTrial(gainTrialsIdx,1); % adjust for picking out data trials later
    
    hitTrial = hitTrial(gainTrialsIdx); % adjust for picking out Hit/Miss later
    
elseif nonGainTrials
    gainTrialsIdx = true(1,size(dFTrials,3));
    gainTrialsIdx(1,trialStart:nthTrial:size(dFTrials,3)) = false;
    dFTrials = dFTrials(:,:,gainTrialsIdx);
    velocityTrials = velocityTrials(:,:,gainTrialsIdx);
    licksTrial = licksTrial(:,:,gainTrialsIdx);
    rwdTrial = rwdTrial(:,:,gainTrialsIdx);
    actlabTrial = actlabTrial(:,:,gainTrialsIdx);
    trialTimes = trialTimes(gainTrialsIdx,1); 
    dataTrial = dataTrial(gainTrialsIdx,1);  % adjust for picking out data trials later
    
    hitTrial = hitTrial(gainTrialsIdx); % adjust for picking out Hit/Miss later
end
% Then pick out either only hits or only misses
if useHitsOnly
    dFTrials = dFTrials(:,:,hitTrial);
    velocityTrials = velocityTrials(:,:,hitTrial);
    licksTrial = licksTrial(:,:,hitTrial);
    rwdTrial = rwdTrial(:,:,hitTrial);
    actlabTrial = actlabTrial(:,:,hitTrial);
    trialTimes = trialTimes(hitTrial);
    dataTrial = dataTrial(hitTrial);
    
    hitTrial = hitTrial(hitTrial); % lastly, reduce hitTrial: all should ==1
elseif useMissOnly
    dFTrials = dFTrials(:,:,~hitTrials);
    velocityTrials = velocityTrials(:,:,~hitTrials);
    licksTrial = licksTrial(:,:,~hitTrials);
    rwdTrial = rwdTrial(:,:,~hitTrials);
    actlabTrial = actlabTrial(:,:,~hitTrials);     
    trialTimes = trialTimes(~hitTrials);
    dataTrial = dataTrial(~hitTrials);
    
    hitTrial = hitTrial(~hitTrials); % lastly, reduce hitTrial: should ==0 since we only collect misses
end

% Mean and SEM over trials
meandF  = nanmean(dFTrials,3);
semdF  = nansem(dFTrials,3);

% Return only trials with data (remove nan trials of other tracks or categories)
if dataTrialsOnly
    dFTrials = dFTrials(:,:,logical(dataTrial));
    velocityTrials = velocityTrials(:,:,logical(dataTrial));
    licksTrial = licksTrial(:,:,logical(dataTrial));
    rwdTrial = rwdTrial(:,:,logical(dataTrial));
    actlabTrial = actlabTrial(:,:,logical(dataTrial));    
    trialTimes = trialTimes(logical(dataTrial),1);
    
    hitTrial = hitTrial(logical(dataTrial));
end

end