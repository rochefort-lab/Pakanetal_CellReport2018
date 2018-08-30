% Turn binned recs into samples to do stats with.
%
% Inputs
% ------
% signal : cell array
%     Input should have `nRec` elements, each of which is an array sized
%     [numRois, numBins]. The signal of interest.
% actlab : cell array
%     Input should have `nRec` elements, each of which is an array sized
%     [1, numBins]. Contains codes indicating which action is performed
%     during each bin.
% ss : cell array
%     Input should have `nRec` elements, each of which is an array sized
%     [1, numBins]. Index of which stimulus in `ustm` is being presented
%     during each bin.
% ustm : cell array
%     Input should have `nRec` elements, each of which is a cell array
%     containing `numStim` elements. The unique stimulus identifiers, as
%     produced by `load_visstim`.
% include_stm : cell array, optional
%     Which ustm values should be included. Values in ustm will be included
%     if the first `n` characters match any element in `include_stm` with
%     length `n`. If empty, all ustm will be included. Default is an empty
%     cell array. To match periods of no stimulus, the string 'none' should
%     be included in `include_stm`. Matches are case-sensitive.
% exclude_stm : cell array
%     Which ustm values should be excluded (even if they match
%     include_stm). If empty, no stim are excluded. Default is {'~'}, which
%     excludes only buffering between stimuli.
%
% Outputs
% -------
% signal : array [numRois, numSamples]
%     The signal of interest.
% actlab : array [1, numSamples]
%     Which action is being performed during each sample, following the
%     code given in the input.
% ss_string : cell array [1, numSamples]
%     Which stimulus is being presented during each sample, in cannonical
%     form.
%
% Notes
% -----
% The stimulus strings in `ustm` are cannonised with `cannonise_ustm` so
% they include the stimulus properties but not the instance number, i.e.
% what is on screen but not which occurance it is.
%
% See also: binByActionAndStim, cannonise_ustm.

function [signal, actlab, ss_string] = binnedrec2samples(...
    signal, actlab, ss, ustm, include_stm, exclude_stm)

% Default inputs ----------------------------------------------------------
if nargin < 5
    include_stm = {};
end
if nargin < 6
    exclude_stm = {'~'};
end

% Input handling ----------------------------------------------------------
assert(iscell(signal), 'Signal should be a cell array.');
assert(iscell(actlab), 'Action labels should be a cell array.');
assert(iscell(ss)    , 'Stimulus labels should be a cell array.');
assert(iscell(ustm)  , 'Stimulus labels should be a cell array.');

nRec = numel(signal);
assert(nRec==numel(actlab), 'Inconsistent number of recs for signal vs actlab.');
assert(nRec==numel(ss)    , 'Inconsistent number of recs for signal vs ss.');
assert(nRec==numel(ustm)  , 'Inconsistent number of recs for signal vs ustm.');


% Main --------------------------------------------------------------------

% Turn our stimulus map into a cell array of strings, in canonical form and
% excluding the instance number
ss_string = cell(size(ss));
for iRec=1:numel(ss)
    ss_string{iRec} = cannon_vstim_strings(ss{iRec}, ustm{iRec});
end

% Concatenate all the samples along dimension 2
signal = cat(2, signal{:});
actlab = cat(2, actlab{:});
ss_string = cat(2, ss_string{:});

% Reduce down to just whitelisted ustm values
if ~isempty(include_stm)
    is_in_whitelist = false(size(ss_string));
    for iWhite = 1:numel(include_stm)
        is_this_entry = strncmp(include_stm{iWhite}, ss_string, ...
            length(include_stm{iWhite}));
        is_in_whitelist = is_in_whitelist | is_this_entry;
    end
    signal = signal(:, is_in_whitelist);
    actlab = actlab(:, is_in_whitelist);
    ss_string = ss_string(:, is_in_whitelist);
end
% Remove blacklisted ustm values
if ~isempty(exclude_stm)
    is_in_blacklist = false(size(ss_string));
    for iBlack = 1:numel(exclude_stm)
        is_this_entry = strncmp(exclude_stm{iBlack}, ss_string, ...
            length(exclude_stm{iBlack}));
        is_in_blacklist = is_in_blacklist | is_this_entry;
    end
    signal = signal(:, ~is_in_blacklist);
    actlab = actlab(:, ~is_in_blacklist);
    ss_string = ss_string(:, ~is_in_blacklist);
end

end


function ss_string = cannon_vstim_strings(ss, ustm)

% Parameters --------------------------------------------------------------
USTM_NO_STIM = 'none';

% Main --------------------------------------------------------------------

% Remove the instance count from unique stimulus IDs
ustm_canon = cannonise_ustm(ustm);

ss_string = cell(size(ss));
% Allocate canonical ustm when there is no stimulus
ss_string(ss==0) = repmat({USTM_NO_STIM}, sum(ss==0), 1);
% Allocate canonical ustm for everything with a stimulus
ss_string(ss~=0) = ustm_canon(ss(ss~=0));

end