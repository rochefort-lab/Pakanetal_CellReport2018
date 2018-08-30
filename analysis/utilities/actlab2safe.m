% Reduces actlab down to only buffered periods of "safe" activity.
%
% Inputs
% ------
% tt : array [1, numTimePoints, numRecs]
%     Time (in seconds) of each sample point.
% actlab : array [1, numTimePoints, numRecs]
%     Which action the animal is performing at each point in time. Encoding
%     should be as per `detect_action`: 0=stationary, 1=twitch, 2=twitch,
%     3=locomoting.
% motion_mode : string, optional
%     Which method to use when cutting motion down to safe. The following
%     options are available:
%       'C' : Only the middle section of locomotion remains. First and last
%           1 second is removed.
%       'CD' : Only the first 1 second of locomotion is removed.
%       'BCD' : Motion is left as-is.
%     Default is 'C'.
%
% Outputs
% -------
% actlab : array [1, numTimePoints, numRecs]
%     Which action the animal is performing at each point in time. -1
%     everywhere except from during the parts of stationary/locomotion
%     epochs deemed "safe". Safe periods of stationary activity will be
%     denoted with `0`, and safe periods of locomotion will be denoted
%     with `3`.
%
% See also motion_onoffset, motionless_onoffset.

function actlab = actlab2safe(tt, actlab, motion_mode)

% Input handling
if nargin<3
    motion_mode = 'C';
end

switch lower(motion_mode)
    case 'c'
        onset_col = 3;
        offset_col = 3;
    case 'cd'
        onset_col = 3;
        offset_col = 4;
    case 'bcd'
        onset_col = 2;
        offset_col = 4;
    otherwise
        error('Unrecognised motion_mode: %s', motion_mode)
end

% Remember the original input
actlab_orig = actlab;

% Make sure input is okay
if any(actlab_orig(:)==-1)
    error('Actlab input should not contain -1 values.');
end

% Initialise a new actlab made of -1 everwhere, to fill in as appropriate
actlab = -1 * ones(size(actlab_orig));

% Fill in the data for motionless epochs
for iRec=1:size(actlab, 3)
    % Find the start/stop indexes
    [still_onsets, still_offsets] = motionless_onoffset(...
        actlab_orig(:, :, iRec), tt(:, :, iRec));
    % Update each epoch
    for iEpoch=1:numel(still_onsets)
        actlab(:, still_onsets(iEpoch):still_offsets(iEpoch), iRec) = 0;
    end
end

% If we are just keeping motion epochs as they are, do that
if strcmpi(motion_mode, 'bcd')
    actlab(actlab_orig==3) = 3;
    return;
end

% Fill in the data for locomotion epochs
for iRec=1:size(actlab, 3)
    % Find the start/stop indexes
    [loco_onsets, loco_offsets] = motion_onoffset(...
        actlab_orig(:, :, iRec), tt(:, :, iRec));

    % Handle one of the column being missing
    onsets = min(loco_onsets(onset_col:offset_col, :), [], 1);
    offsets = max(loco_offsets(onset_col:offset_col, :), [], 1);

    % Update each epoch
    for iEpoch=1:size(loco_onsets,2)
        if isnan(onsets(iEpoch)) || isnan(offsets(iEpoch))
            continue;
        end
        actlab(1, onsets(iEpoch):offsets(iEpoch), iRec) = 3;
    end
end

end