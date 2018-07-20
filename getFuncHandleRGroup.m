% Look up the handle for a function which works on a set of recordings.
%
% The function handle is taken from a subfunction within this .m file.
% All functions returned have the format
%     X = func(signal, tt, auxlab, ped, vel, acc)
% for stimulus-independent or
%     [X, ustm] = func(signal, tt, auxlab, ped, vel, acc, stimTimes, ustm)
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
%
% Created by: Rochefort Lab
%
%
% ======================================================================= %
%                           LIST OF FUNCTIONS                             %
% ======================================================================= %
% (Ctrl+D to jump to function in this file)
%__________________________________________________________________________
%   OSI
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
%                       ANALYSIS FUNCTIONS                                %
% ======================================================================= %
% Description of possible inputs
%__________________________________________________________________________
% Inputs
% ------
%   signal : array [nROIs, nSamples, nRecs]
%       Fluorescence activity traces per ROI for each trial and track
%       nSamples = tot number of samples during the recording
%   tt : array [1, nSamples, nRecs]
%       Time of samples since the start of the recording (seconds)
%   auxlab : array [1, nSamples, nRecs]
%     Action labels of time points, in the format given by detect_action
%   stimTimes: cell array [nRecs, 1]
%       Each cell contains an array [1, nStim]
%   ustm: cell array [nRecs, 1]
%       Each cell contains a cell array [1, nStim]
%__________________________________________________________________________





function [X, stimFlag] = OSI(signal, tt, ~, ~, ~, ~, stimTimes, ustm)
% Calculate OSI with pref-orth/pref+orth method. 
% Do not separate out actions. 
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 11]   For each ROI (Columns):
%       1  - OSI magnitude of the average rec
%       2  - OSI pref angle of the average rec
%       3  - Average OSI magnitude over recs 
%       4  - SEM OSI magnitude over recs 
%       11 - Std OSI magnitude over recs 
%       5  - Mode(or mean for CV) of preferred orientation over recs
%       6  - Std of preferred orientation over recs
%       7  - Average OSI magnitude over recs with preferred ori (sucessful response)
%       8  - SEM OSI magnitude over recs with preferred ori (sucessful response)
%       9  - Percentage of trials with sucessful response (columns 7&8)
%       9  - Total number of trials
%
% see also: calcOSI
%__________________________________________________________________________

stimFlag = [];

% 1) Convert function input to relevant Data 

[ dataDf_stimMean, ustm, stimAngle ] = prepdataforfunction_RGroup('OSI', signal, tt, stimTimes, ustm);
% - Here: extracted data for Grating stim           

% 2) Analysis 
%
% Format of Data for analysis
% ===========================
%   dataDf_stimMean : cell array [nRecs, 1] 
%       Each cell contains the average signal of each ROI within a period
%       of stimulation for each rec 
%   ustm : cell array [nRecs, 1]
%   stimAngle : cell array [nRecs, 1] Stim angle protocol for each rec

nRecs = length(dataDf_stimMean);
nROIs = size(dataDf_stimMean{1},1);

% calculate OSI with all data
stimMeanALL = reshape(cat(2,dataDf_stimMean{:}), nROIs, []);
[prefOriALL, prefSIALL] = calcOSI(stimMeanALL, [ustm{:}], [stimAngle{:}], true, false);

% calculate OSI for each rec
prefOri = nan(nROIs, nRecs);
prefSI = nan(nROIs, nRecs);
for iRec = 1:nRecs
    % calculate OSI per Rec
    [prefOri(:,iRec), prefSI(:,iRec)] = calcOSI(dataDf_stimMean{iRec}, ustm{iRec}, stimAngle{iRec}, true, false);
end


% -- Output ---------------------------------------------------------------
% Return total mean OSI magnitude for all recs meaned first (column 1) and OSI pref angle for all recs meaned first (column 2)
X(:,1) = prefSIALL;
X(:,2) = prefOriALL;
% Return total mean OSI magnitude over recs, meaned after (column 3) and sem (column 4).
X(:,3) = nanmean(prefSI,2);
X(:,4) = nansem(prefSI,2);
X(:,11) = nanstd(prefSI,[],2);
% Return mode(or mean for CV) (column 5) and std (column 6) OSI preferred orientation across recs
X(:,5) = mode(prefOri,2);
% X(:,5) = nanmean(prefOri,2); % for circular variance get mean across trials
X(:,6) = nanstd(prefOri,0,2);
% Return total mean OSI magnitude over recs with preferred ori (sucessful response), meaned after (column 7) and sem (column 8).
for iROI = 1:nROIs
    is_prefOriRec = prefOri(iROI,:)==prefOriALL(iROI);
%     prefOriRec = prefOri(iROI,:)<=(prefOriALL(iROI)+5) & prefOri(iROI,:)>=(prefOriALL(iROI)-5); % for circular variance, must be within 5 degrees of preferred ori across all trials
    X(iROI,7) = nanmean(prefSI(iROI,is_prefOriRec),2);
    X(iROI,8) = nansem(prefSI(iROI,is_prefOriRec),2);
    % Return the in column 9 the percent of trials used for calculations in column 7 and 8 (i.e. % of trials with preferred orientation), and in column 10 the total number of trials
    X(iROI,9) = sum(double(is_prefOriRec))/nRecs; % percent trials used
end
X(:,10) = nRecs; % total number of trials

end

