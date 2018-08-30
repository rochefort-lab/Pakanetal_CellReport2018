% Bins recordinngs so they are categorised by action and stimulus
%
% Inputs
% ------
% tt : array [1, numTimePoints, numRecs]
%     Time at which samples are measured, in seconds.
% actlab : array [1, numTimePoints, numRecs]
%     Action labels of time points, in the format given by detect_action.
% stimTimes : cell array
%     Time of stimulus onsets. Should contain `numRecs` many cells, each
%     sized [1, numStim]. This should not be a 3-D array, because the
%     number of of stimuli can be different for each recording.
% signal : array [numRois, numTimePoints, numRecs]
%     The signal of interest.
% mode : string, optional
%     Controls what method is used to determine the bin output. Should be
%     one of the following:
%       'filter' : `signal` is downsampled to a lower sampling rate with
%           interval equal to `binDur` (i.e. a sampling frequency of
%           `1/binDur` Hz), at which point samples are assumed to be
%           independent from each other.
%       'bin-mean' : `signal` is broken into bins with duration `binDur`
%           and at least `binIntv` spacing between them, then the average
%           of `signal` during these bins is the output.
%       'bin-change' : `signal` is broken into bins with duration `binDur`
%           and at least `binIntv` spacing between them, then the output is
%           the difference between the average of `signal` in the last
%           0.1s of the bin and the first 0.1s.
%     Default is 'filter'.
% binDur : double, optional
%     Duration of sample bins, in seconds. Default is 1.75.
% binIntv : double, optional
%     Minimum inter-bin interval, in seconds. Default is 0.
% maxOneBinPerStim : bool, optional
%     Whether to allow multiple bins to come from the same stimulus if it
%     is longer than `binDur`. NB: multiple bins can always come from
%     periods without any stimulation. Default is `true`.
%
% Outputs
% -------
% tt : cell array [1, 1, numRecs]
%     Each cell is an array, sized [1, numBins]. Time at which each bins
%     begins, in seconds.
% actlab : cell array [1, 1, numRecs]
%     Each cell is an array, sized [1, numBins]. Action label for each bin.
% ss : cell array [numRecs]
%     Each cell is an array, sized [1, numBins]. Index of the stimulus
%     which is present during the bin, pertaining to `stimTimes`. An index
%     of 0 in `ss` indicates the bin is before the first stimulus.
% signal : cell array [1, 1, numRecs]
%     Each cell is an array, sized [numRois, numBins]. Bin values. Either
%     downsampled version of input, or bin mean, or difference between
%     start and end of bin, depending on `mode`.
%
% See also: binnedrec2samples, getFuncHandleRGroup>AnovaVisuallyResponsive.

function [tt, actlab, ss, signal] = binByActionAndStim(...
    tt, actlab, stimTimes, signal, mode, binDur, binIntv, maxOneBinPerStim)

% We want to find each the edge of each region.
% Regions are the same action and same stimulus throughout.
% Then find the maximum nubmer of bins with duration binDur which fit in
% the region.
% The bins are arranged so they are spaced apart as far as possible.
% The varargin variables are averaged within each bin.
% The bined variables are returned, along with an array indicating which
% action and stimulus they are within.

% Alternatively, we downsample the data to 1/binDur rate, and then
% identify which of the resulting datapoints are in each region.

% Settings ----------------------------------------------------------------
stimLag   = 0;         % Expected lag from stimulus to response

% Default inputs ----------------------------------------------------------
if nargin<5
    mode = 'filter'; % which method to use
end
if nargin<6
    binDur = 1.8;  % Bins this long (seconds)
end
if nargin<7
    binIntv = 0;  % Bins with this spacing interval between them
end
if nargin<8
    maxOneBinPerStim = true;
end

% Input handling ----------------------------------------------------------
assert(iscell(stimTimes), 'Need cell input for stimTimes')
% Ensure inputs are cells
if ~iscell(tt)
    tt     = num2cell(tt, [1,2]);
end
if ~iscell(actlab)
    actlab = num2cell(actlab, [1,2]);
end
if ~iscell(signal)
    signal = num2cell(signal, [1,2]);
end
nRec = numel(stimTimes);
assert(nRec==numel(tt)    , 'Number of recordings are inconsistent');
assert(nRec==numel(actlab), 'Number of recordings are inconsistent');
assert(nRec==numel(signal), 'Number of recordings are inconsistent');


% Main --------------------------------------------------------------------
if strcmp(mode, 'filter')
    % For each tt, which stim is it in?
    ss = cell(size(stimTimes));
    for iRec=1:nRec
        [signal{iRec}, tt{iRec}, actlab{iRec}] = decimate_rec(1/binDur, ...
            signal{iRec}, tt{iRec}, actlab{iRec});
        ss{iRec} = stimTimes2stimmap(stimTimes{iRec}, tt{iRec});
    end
    return;
    
elseif strncmp(mode, 'bin-', 4)
    
    ss = cell(size(stimTimes));
    
    for iRec=1:nRec
        [tt{iRec}, actlab{iRec}, ss{iRec}, signal{iRec}] = bin_single_rec(...
            tt{iRec}, actlab{iRec}, stimTimes{iRec}, signal{iRec}, ...
            binDur, binIntv, mode(5:end), maxOneBinPerStim);
        
    end
    
else
    error('Unrecongised mode: %s', mode);
end

end


function [tt_out, actlab_out, ss_out, signal_out] = bin_single_rec(...
    tt, actlab, stimTimes, signal, bin_dur_sec, bin_intv_sec, method, ...
    maxOneBinPerStim)
% Bin a single recording

% Parameters --------------------------------------------------------------
change_win_sec = 0.1;


% Default inputs ----------------------------------------------------------
if nargin<6 || isempty(bin_intv_sec)
    bin_intv_sec = 0;
end
if nargin<7 || isempty(method)
    method = 'mean';
end
if nargin<8
    maxOneBinPerStim = true;
end


% Input handling ----------------------------------------------------------
assert(~iscell(tt)       , 'Need array input for tt');
assert(~iscell(actlab)   , 'Need array input for actlab');
assert(~iscell(stimTimes), 'Need array input for stimTimes');
assert(~iscell(signal)   , 'Need array input for signal');

% Ensure inputs are cells
assert(ndims(tt)    <3   , 'Single rec input for tt');
assert(ndims(actlab)<3   , 'Single rec input for actlab');
assert(ndims(signal)<3   , 'Single rec input for signal');
assert(length(stimTimes)==numel(stimTimes), 'StimTimes must be a vector');

tt_intv = median(diff(tt));
bin_dur_idx    = round(bin_dur_sec ./ tt_intv);
bin_intv_idx   = round(bin_intv_sec ./ tt_intv);
change_win_idx = round(change_win_sec ./ tt_intv);

if strcmpi(method, 'change')
    assert(bin_dur_idx > 2*change_win_idx, 'Bin is too short');
end

if maxOneBinPerStim
    stmToMultiBin = 0;
else
    stmToMultiBin = 0:numel(stimTimes);
end


% Main --------------------------------------------------------------------

% Initialise outputs
tt_out     = nan(size(tt,1), 0);
actlab_out = nan(size(actlab,1), 0);
ss_out     = nan(size(ss,1), 0);
signal_out = nan(size(signal,1), 0);

% Pad stim times so we know the start of the pre-stim period, and
% the end of the final stimulus (which is the end of the rec)
paddedStimTimes = [tt(1); stimTimes(:); tt(end)];


% Prep for traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ED edit start
if size(signal,2)<size(tt,2);
    R=size(tt,2);
    rem=size(signal,2);
% Check if signal is divisible by factor
        padNum = R-rem;
        if padNum>1/mean(diff(tt));
            error('wheel and signal differ by more than 1 second (%d frames)',padNum)
        end
        warning('Uneven number of frames in recording, padding with last frame value for %d frames',padNum)

    for i=1:padNum
        signal = [signal, signal(:,end,:)]; %#ok<AGROW>
    end
end
% Prep for traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ED edit end


% Bin each rec
for iStm=0:numel(stimTimes)
    
    isThisStim = tt >= paddedStimTimes(iStm+1) & ...
                 tt <= paddedStimTimes(iStm+2);

    tt_stm = tt(:, isThisStim);
    al_stm = actlab(:, isThisStim);
    sg_stm = signal(:, isThisStim);
    
    % Find where action changes
    idx_change_al_stm = find(diff(al_stm));
    % We pad with 0 and length of the vector so have reference points for
    % the start of the first and end of the final action.
    idx_change_al_stm = [0, idx_change_al_stm, length(al_stm)];
    
    nAct = length(idx_change_al_stm)-1;
    % If we are only doing one bin, that bin is at the start of the
    % stimulus, even if there are multiple actions present. That is to say,
    % we will only use the action which occurs first.
    if ~ismember(iStm, stmToMultiBin)
        nAct = min(1, nAct);
    end
    
    for iAct=1:nAct
        
        % Check how long this action is present for 
        al_dur_idx = idx_change_al_stm(iAct+1) - idx_change_al_stm(iAct);
        
        % We need to have at least the bin length and the bin interval
        % between each bin
        nBin = floor(al_dur_idx / (bin_dur_idx+bin_intv_idx));
        
        % If there isn't enough duration, check if we can still fit one in
        % without the spacing present
        if nBin==0 && al_dur_idx >= bin_dur_idx
            nBin = 1;
        end
        
        if nBin==0;
            continue;
        end
        
        if ~ismember(iStm, stmToMultiBin)
            nBin = min(1, nBin);
            bin_start_offsets = 1;
        else
            spare_idx = al_dur_idx - nBin*bin_dur_idx;
            spare_idx_per_bin = spare_idx/nBin;
            bin_start_offsets = (spare_idx_per_bin/2 + 1) : ...
                (bin_dur_idx + spare_idx_per_bin) : al_dur_idx;
            bin_start_offsets = round(bin_start_offsets);
        end
        
        bin_starts = idx_change_al_stm(iAct) + bin_start_offsets;
        
        for iBin=1:nBin
            ii = bin_starts(iBin) + (0:(bin_dur_idx-1));
            tt_out(:,end+1)     = tt_stm(:,bin_starts(iBin));
            actlab_out(:,end+1) = al_stm(:,floor(bin_starts(iBin)+bin_dur_idx/2));
            ss_out(:,end+1)     = iStm;
            if strcmpi(method, 'mean')
                signal_out(:,end+1) = mean(sg_stm(:,ii),2);
            elseif strcmpi(method, 'change')
                ii_on  = bin_starts(iBin) + (0:(change_win_idx-1));
                ii_off = bin_starts(iBin) - (0:(change_win_idx-1)) + bin_dur_idx-1;
                signal_out(:,end+1) = mean(sg_stm(:,ii_off),2) ...
                                      - mean(sg_stm(:,ii_on),2);
            else
                error('Unrecognised method: %s', method);
            end 
        end
    end
end

end