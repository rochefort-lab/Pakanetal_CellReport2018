function [signal, actlab, stm, correlation] = generic_sampler(...
    signal, tt, actlab, stimTimes, ustm, include_stm, exclude_stm)

% Default inputs
if nargin<6
    include_stm = {};
end
if nargin<7
    exclude_stm = {'~'};
end

% Parameters
motion_mode = 'BCD';  % only drop "un-safe" stationary periods, leave loco
binDur = 1.0;  % bin duration, in seconds (no longer than shortest stimulus period)
binIntv = 2.0;  % bin interval, in seconds 
maxOneBinPerStim = false;  % allow multiple samples to come from one epoch

% Cut actions down to safe actions
actlab = actlab2safe(tt, actlab, motion_mode);
% Bin the data
[tt, actlab, ss, signal] = binByActionAndStim(...
    tt, actlab, stimTimes, signal, 'bin-mean', binDur, binIntv, ...
    maxOneBinPerStim);

% Order our binned samples nicely
[signal, actlab, stm] = binnedrec2samples(...
    signal, actlab, ss, ustm, include_stm, exclude_stm);

% Remove uninteresting actions
is_nice_action = ismember(actlab, [0, 3]);

% Remove unwanted samples
signal = signal(:, is_nice_action);
actlab = actlab(:, is_nice_action);
stm = stm(:, is_nice_action);

if nargin >= 4
    correlation = binsettings2correlation(binDur, binIntv);
end

end