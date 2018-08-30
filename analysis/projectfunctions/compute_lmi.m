% Computes an unbiased estimate of the locomotion modulation index
%
% Computes the LMI by breaking the signal down into (roughly) independent
% samples, which each occur during a single visual stimulus and a single
% activity of the animal. These bins have a duration of 1 second and at
% least 1 second in between them. Stationary epochs are restricted to be at
% least 3 seconds after the last locomotion epoch (and end 0.2 seconds
% before the next one; see motionless_onoffset). The LMI is computed for
% each unique stimulus (stationary grating orientation 0 degrees, moving
% grating 90 degrees, etc) by taking the relative difference between
% average signal during this stimulus under locomotion and stationary
% activity. From this, the average LMI for each stimulus category
% (stationary, moving) is computed. These are then averaged to find an
% unbiased estimate of the LMI under all conditions. Note that uniform
% stimuli (which have ustm beginning with 'U') are ignored.
%
% Inputs
% ------
% signal : array [numRois, numTimePoints, numRecs]
%     The signal which may or may not be modulated by locomotion.
% tt : array [1, numTimePoints, numRecs]
%     Time at which each sample occurs (in seconds).
% actlab : array [1, numTimePoints, numRecs]
%     Encoding of which action is being performed at each point in time.
% stimTimes : cell array [numRecs]
%     Time at which each stimulus begins (in seconds). Each cell should
%     contain an array of doubles, with an element per stimulus.
% ustm : cell array [numRecs]
%     Unique stimulus identifiers for each stimulus whose timing is
%     recorded in `stimTimes`. Each cell should contain a cell array of
%     strings, with an element per stimulus.
% num_bootstraps : int, optional
%     How many bootstraps to perform. Default is 1000. Note: No bootstraps
%     will be calculated if the number of outputs is less than 2.
% bootstrap_prop : double, optional
%     Proportion of the data to include in each bootstrap subset of the
%     data samples. If empty, `binsettings2correlation` is used to
%     determine how many samples should be included in the bootstraps.
%     Default is []. Note that bootstrap resamples are drawn without
%     replacement.
% include_stm : cell array, optional
%     Which ustm values should be included. Values in ustm will be included
%     if the first `n` characters match any element in `include_stm` with
%     length `n`. If empty, all ustm will be included. Default is an empty
%     cell array. To match periods of no stimulus, the string 'none' should
%     be included in `include_stm`. Matches are case-sensitive.
% exclude_stm : cell array
%     Which ustm values should be excluded (even if they match
%     include_stm). If empty, no stim are excluded. Default is an empty
%     cell array.
% rng_seed : int, optional
%     The seed to initialise the random number generator with before
%     bootstrapping. Default is `570`.
%
% Outputs
% -------
% lmi : array [numRois, 1]
%     Unbiased estimate of LMI of `signal`, conditioned on the stimuli
%     `ustm` which occur at `stimTimes`.
% lmi_btsp : array [numRois, numBootstraps]
%     Bootstrap outputs corresponding to `lmi` when only `bootstrap_prop`
%     of the samples are included.
%
% See also actlab2safe, binByActionAndStim, binnedrec2samples, samples2lmi.

function [lmi, lmi_btsp] = compute_lmi(...
    signal, tt, actlab, stimTimes, ustm, num_bootstraps, bootstrap_prop, ...
    include_stm, exclude_stm, rng_seed)
% Computes an unbiased estimate of the locomotion modulation index

% Default inputs
if nargin<6
    num_bootstraps = 1000;
end
if nargin<7
    bootstrap_prop = [];
end
if nargin<8
    include_stm = {};
end
if nargin<9
    exclude_stm = {};
end
if nargin<10
    rng_seed = 570;
end

nRoi = size(signal, 1);

[signal, actlab, stm, cor] = generic_sampler(...
    signal, tt, actlab, stimTimes, ustm, include_stm, exclude_stm);

if isempty(bootstrap_prop) && num_bootstraps > 0
    % Use the geometric series to work out how much each datapoint gets counted
    %          S = 1   + x   + x^2 + x^3 + ...
    %        S*x = x   + x^2 + x^3 + ...
    %    S - S*x = 1
    %   S(1 - x) = 1
    %          S = 1 / (1 - x)

    % We need to correct for this, so we are interested in using 1/S as the
    % proportion
    %   prop = 1 / S = 1 - x
    bootstrap_prop = 1 - cor;
end

% Deal with no samples available
if isempty(signal)
    lmi = nan(nRoi, 1);
    lmi_btsp = nan(nRoi, num_bootstraps);
    return;
end

% Compute an unbiased estimate of the lmi from our samples
lmi = samples2lmi(signal, actlab, stm);

if nargout<2
    % Don't need to bootstrap, so stop now
    return;
end

% Set the RNG, for repeatability
prev_state = rng(rng_seed, 'twister');
% Initialise
lmi_btsp = nan(size(signal,1), num_bootstraps);
% Do the bootstrapping
for iBtsp = 1:num_bootstraps
    % Draw random resampling, with replacement
    indices = randi(size(signal,2), floor(size(signal,2) * bootstrap_prop), 1);
    % Compute LMI from this resampling
    lmi_btsp(:, iBtsp) = samples2lmi(...
        signal(:, indices), actlab(:, indices), stm(:, indices));
end

% Reset the RNG
rng(prev_state);

end