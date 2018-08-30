% Organize fluorescence dF/F0 signals per track for each ROI.
% 
% Each track has a set of trials. For each ROI calculate dF/F0 for all 
% trials of each track. Signals are binned either by distance or time steps
% 
% Inputs
% ------
%   signal : array [nROIs, nSamples]
%       Fluorescence activity traces per ROI for each trial and track
%       nSamples = tot number of samples during the recording
%   tt : vector [1, nSamples]
%       Time of samples since the start of the recording (seconds)
%   VRlab : structure with fields:
%       Each field is a vector of size [1, nSamples]
%       - trackNum
%       - reward
%       - trialNum
%       - licks
%       - actlab
%   displmnt : vector [1, nSamples]
%       Cumulative total distance travelled (displacement) since the start 
%       of the recording (cm)
%   vel : vector [1, nSamples]
%       Instantaneous velocity (cm/s)
%   lagTime : vector [1, nSamples]
%       Instantaneous frame by frame lag in VR (seconds)
%        
% Inputs (optional)
% ----------------
%   binnr :  vector [1, numTracks]   (Default: [48, 48, 64, 64, 96, 96])
%       Each element is the number of bins for each track number(see below) 
%       Default for tracks of length 120, 160 or 240 cm
%       Bins are 2.5cm
%   trackNums : vector [1, numTracks]   (Default: [1,2,3,5,7,8])
%       Requested tracks (VR world number, set by Virmen)
%   binBy : string   (Default: 'distance')   
%       Possible values = 'distance' , 'time'
%   dataTrialsOnly : boolean   (Default: false)
%       Return only trials with data (remove nan trials of other tracks or
%       categories)
%   useHitsOnly : boolean   (Default: false)
%       Some descript
%   useMissOnly : boolean   (Default: false)
%       Some descript
%
% Outputs
% -------
%   sortROIIndx :
%       Description
%   dFTrials :
%       Description
%   normdF :
%       Description
%   meandF :
%       Description
%   semdF :
%       Description
%   sortdF :
%       Description
%   hitTrial : cell [numTracks, 1]
%       Declares for each trial of each track if the animal got a reward.
%       Each cell contains a boolean vector [numTotTrials_track, 1] where
%           numTotTrials_track = Total number of trials for the track
%   velTrials : cell [numTracks, 1]
%       Declares for each trial of each track the animals velocity along the track.
%       Each cell contains a vector [numTotTrials_track, nBins] where
%           numTotTrials_track = Total number of trials for the track and
%           nBins = number of bins along length of the track (e.g. binned by 2.5 cm)
%   trialTimes : cell [numTracks, 1]
%       Declares for each trial of each track the total trial length in seconds.
%       Each cell contains a vector [numTotTrials_track, 1] where
%           numTotTrials_track = Total number of trials for the track
%   rwdTrials : cell [numTracks, 1]
%       Declares for each trial of each track where the animal got a reward (1 or 2) along the track (binned).
%       Each cell contains a vector [numTotTrials_track, nBins] where
%           numTotTrials_track = Total number of trials for the track and
%           nBins = number of bins along length of the track (e.g. binned by 2.5 cm)
%   licksTrials : cell [numTracks, 1]
%       Declares for each trial of each track if the animal licked.
%       Each cell contains a vector [numTotTrials_track, nBins] where
%           numTotTrials_track = Total number of trials for the track and
%           nBins = number of bins along length of the track (e.g. binned by 2.5 cm)
%   actlabTrials : cell [numTracks, 1]
%       Declares for each trial of each track the action label of the animal.
%       Each cell contains a vector [numTotTrials_track, nBins] where
%           numTotTrials_track = Total number of trials for the track and
%           nBins = number of bins along length of the track (e.g. binned by 2.5 cm)
%
% See also: select_tracks, bin_dF_distance

function [sortROIIndx, dFTrials, normdF, meandF, semdF, sortdF, hitTrial,...
    velTrials, trialTimes, rwdTrials, licksTrials, actlabTrials] =...
    VRdata2tracks(signal, tt, VRlab, displmnt, vel, lagTime, ...
        binnr, trackNums, binBy, dataTrialsOnly, useHitsOnly, useMissOnly)

if nargin<12 || isempty(useMissOnly)
    useMissOnly = false;
end
if nargin<11 || isempty(useHitsOnly)
    useHitsOnly = false;
end
if nargin<10 || isempty(dataTrialsOnly)
    dataTrialsOnly = false;
end
if nargin<9 || isempty(binBy)
    binBy = 'distance';
end
if nargin<8 || isempty(trackNums)
    trackNums = [1,2,3,5,7,8];
end
if nargin<7 || isempty(binnr)
    binnr = [48, 48, 64, 64, 96, 96];
end

numTracks = length(trackNums);
nROIs = size(signal,1);

% reassemble into original behaviour_aligned organization
VRdata = [tt; displmnt; lagTime; vel; cell2mat(struct2cell(VRlab))];

% Initialize variables for each track
VRdataTracks = cell(numTracks,1);
dFdataTracks = cell(numTracks,1);
dFTrials = cell(numTracks,1);
meandF = cell(numTracks,1);
semdF = cell(numTracks,1);
hitTrial = cell(numTracks,1);
normdF = cell(numTracks,1);
sortROIIndx = nan(nROIs,numTracks);
sortdF = cell(numTracks,1);
trialTimes = cell(numTracks,1); % Get total trial length (in sec) for each trial, separated by track ID (*Note, this is {1 x trial}, not binned).
velTrials = cell(numTracks,1);  % Get velocity (in cm/sec) binned across trach length and for each trial, separated into cells by track ID number.
rwdTrials = cell(numTracks,1);
licksTrials = cell(numTracks,1);
actlabTrials = cell(numTracks,1);

for iTrack = 1:length(trackNums)
    
    [VRdataTracks{iTrack}, dFdataTracks{iTrack}] = select_tracks(VRdata, 'all', trackNums(iTrack), signal);
    
    if ~isempty(VRdataTracks{iTrack})
        [dFTrials{iTrack}, meandF{iTrack}, semdF{iTrack}, hitTrial{iTrack},...
            velTrials{iTrack}, trialTimes{iTrack}, rwdTrials{iTrack}, licksTrials{iTrack}, actlabTrials{iTrack}]...
        = bin_dF_distance(VRdataTracks{iTrack}, dFdataTracks{iTrack}, binnr(iTrack), binBy, dataTrialsOnly, useHitsOnly, [], [], useMissOnly);
        
        normdF{iTrack} = meandF{iTrack}./repmat(max(meandF{iTrack},[],2),1,size(meandF{iTrack},2));
        
        [~, maxXIndx] = max(normdF{iTrack}, [], 2);
        [~, sortROIIndx(:,iTrack)] = sort(maxXIndx);
        dataTable = normdF{iTrack};
        sortdF{iTrack} = dataTable(sortROIIndx(:,iTrack), :);
    end
end

end