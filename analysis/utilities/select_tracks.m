% Select data for given tracks.
%
% Return the data of all trials corresponding to the requested tracks.
% Optional: also returns the correponding fluorescence data
% 
% Inputs
% ------
%   VRdata : array [9, nSamples] 
%       nSamples = tot number of samples during the recording
%       Rows in VRdata:
%       1 - Time of samples since the start of the recording (seconds)
%       2 - Distance travelled since the start of the recording (cm)
%       3 - Instantaneous frame by frame lag in VR (s)
%       4 - Instantaneous velocity (cm/s)
%       5 - trackNum (VR world number - output from Virmen)
%       6 - reward (event detection of reward;  0=none, 1=early, 2=late)
%       7 - trialNum (running trial counter)
%       8 - licks (event detection from lick sensor)
%       9 - actlab (action for animal state - 0=stationary, 1|2=positioning, 3=locomotion)
%   dset_type : string 
%       Possible values = 'all' , 'vel' , 'licks' , 'reward'
%       If none of the above, returns an error
%   track : vector with integers [1, numTrackSelect]
%       Requested track(s) to select among 'trackType'(5th row in 'VRdata')
%       
% Inputs (optional)
% ----------------
%   dfdata :  array [nROIs, nSamples]   (Default: [])
%       If provided, also return deltaF data 
%
% Outputs
% -------
%   VRdataTracks : array [9, nSamplesTracks]
%       Same as input VRdata but only with the frames corresponding to the
%       requested tracks
%   dfdata : array [nROIs, nSamplesTracks]
%               or    return [] if input 'dfdata' is empty or not provided
%
% 
% See also: VRdata2tracks

function [VRdataTracks, dfdata] = select_tracks(VRdata, dset_type, track, dfdata)

doDF = true; % if provided, also return deltaF data
if nargin<4 || isempty(dfdata)
    doDF = false;
end
    
trackNum = VRdata(5,:); % Isolate the track time from the VR data (behaviour_aligned in the hd5 file)
% decide which behaviour data parameters to return (rows)
if strcmp(dset_type,'all')
    % Do nothing
elseif strcmp(dset_type,'vel')
    VRdata = VRdata(4,:);
elseif strcmp(dset_type, 'licks')
    VRdata = VRdata(8,:);
elseif strcmp(dset_type,'reward')
    VRdata = VRdata(6,:);
else
    error('Unrecognised dataset type. Please use either ''all'', ''vel'', ''licks'', or ''reward''');
end

% remove data depending on which track(s) were selected
if length(track)==1 && track ~= -1 
    VRdataTracks = VRdata(:,trackNum==track);
    if doDF
        dfdata = dfdata(:,trackNum==track);
    end
else
    VRdataTracks = VRdata(:,ismember(trackNum, track));
    if doDF
        dfdata = dfdata(:,ismember(trackNum, track));
    end
end

end