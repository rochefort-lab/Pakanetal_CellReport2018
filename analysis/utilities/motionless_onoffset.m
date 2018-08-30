% Finds the indices at which stationary periods begin and end, whilst
% avoiding the end of any non-stationary periods with a buffer of 3 secs.
% 
% Inputs
% ------
% actlab : array [1, numTimePoints]
%     Action labels of time points, in the format given by detect_action.
% tt : array [1, numTimePoints]
%     Times of sampling points (in seconds).
%
% Outputs
% -------
% onsets : array [1, numMotionlessPeriods]
%     Indices of onsets of stationary periods.
% offsets : array [1, numMotionlessPeriods]
%     Indices of offsets of stationary periods.

function [onsets, offsets] = motionless_onoffset(actlab, tt)

% -------------------------------------------------------------------------
% Parameters
pre_avoid_dur    = 0.2; % seconds
post_avoid_dur   = 3.0; % seconds

% -------------------------------------------------------------------------
% Sampling frequency is inverse of the difference in sampling times
fs = 1./mean(diff(tt));

% Convert buffer time durations into indices
pre_avoid_idx  = round(pre_avoid_dur  * fs);
post_avoid_idx = round(post_avoid_dur * fs);

% Need to pad the action labels to find edges at the start/end of recording
pad_al = [NaN, actlab, NaN];

% We want to find stationary
isstill = pad_al==0;

% Find start and end of stationary periods
onsets  = 1 + find(diff(isstill)== 1); % we add 1 so relative to pad_al start
offsets =     find(diff(isstill)==-1);

% Pad all the starts so they are delayed by the appropriate duration
% Except if it is stationary from the very start
onsets(onsets>2) = onsets(onsets>2) + post_avoid_idx;
% Pad all the ends so they avoid the beginning of the next action
% Except if it is stationary to the very end
offsets(offsets<length(pad_al)-1) = offsets(offsets<length(pad_al)-1) - pre_avoid_idx;

% Discard any stationary periods which end before they begin
isGood  = onsets < offsets;
onsets  = onsets(isGood);
offsets = offsets(isGood);

% Take off the padding
onsets  = onsets -1;
offsets = offsets-1;

end
