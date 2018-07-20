% Calculate orientation selectivity index.
%
% Given imaging data for a selection of stimuli, computes the OS
% (or, optionally, DSI - direction selectivity index) for the
% moving grating in the presented stimulus set.
%
% Inputs
% ------
% stimData : double array sized [numRois, numStimuli, ...]
%     Recording data.
% ustm : cell array of strings sized [1, numStimuli]
%     Unique stimulus identities
% stimAngle : vector
%     Angle of stimuli in degrees.
% FLAG_OSI : boolean, optional
%     Whether to compute the OSI or the DSI. If true, OSI is computed.
%     otherwise the DSI is computed instead. Deafult is true.
% useCVmethod : boolean, optional
%     Whether to use circular variance (Ringach et al., 2002, J Neurosci)
%     to calculate OSI. Defaults to true. If set as false then OSI is
%     calculated as (Pref-Orth)/(Pref+Orth).
%
% Outputs
% -------
% prefOri : double array sized [numRois, 1, ...]
%     Prefered orientation of each ROI, in degrees.
% prefSI : double array sized [numRois, 1, ...]
%     Magnitude of the orientation selectivity (OSI or DSI).

function [prefOri, prefSI] = calcOSI(stimData, ustm, stimAngle, FLAG_OSI, useCVmethod)

% Default inputs ----------------------------------------------------------
if nargin<5 || isempty(useCVmethod) 
    useCVmethod = true; % use circular variance to calculate OSI
end
if nargin<4 || isempty(FLAG_OSI)
    FLAG_OSI = true; % for 180degree
end

% Input handling ----------------------------------------------------------
assert(ismatrix(stimAngle), 'stimAngle cant have more than 2 dimensions');
if size(stimAngle,1)>1
    stimAngle = stimAngle';
end
assert(size(stimAngle,1)==1, 'stimAngle must be a vector not a matrix');

assert(size(stimData,2)==numel(stimAngle), 'Input sizes do not match');

% Main --------------------------------------------------------------------

% Check which stimuli are moving (either F for forward or B for backward as
% the first letter in ustm). We assume everything else is static.
isMoving = strncmp(ustm,'F',1) | strncmp(ustm,'B',1);

% Cut down to only the moving stimuli
stimAngle = stimAngle(isMoving);

% Cut down to only the moving stimuli
pp = repmat({':'},1,ndims(stimData));
pp{2} = isMoving;
stimData = stimData(pp{:});


% If we are doing orientation selectivity, double the angle to wrap 180
% degrees into a full circle. Now 90 degrees and 270 degrees are the same
% orientation instead of different directions.
if FLAG_OSI
    stimAngle = 2 * stimAngle;
end

% Use modulo to get principle versions of all angles
stimAngle = mod(stimAngle, 360);
% Check which angles exist
[unique_stimAngle, ~, ic] = unique(stimAngle);

% Give an error if unique_stimAngle is not evenly distributed
da = diff([unique_stimAngle(:); unique_stimAngle(1)+360]);
assert(all(abs(da - mean(da)) < 10^-9), 'Not all angles were sampled');

% Inititalise matrix for average response to each angle
siz = size(stimData);
siz(2) = numel(unique_stimAngle);
average_responses = nan(siz);
% Take the average for each angle
for iAngle=1:numel(unique_stimAngle)
    % Need to do slicing of unknown number of dimensions on left and right
    % hand sides of equation
    lft = repmat({':'},1,ndims(average_responses));
    lft{2} = iAngle;
    rgt = repmat({':'},1,ndims(stimData));
    rgt{2} = ic==iAngle;
    % Take the avearage over all occurances of this angle and put it in
    % the holding array
    average_responses(lft{:}) = nanmean(stimData(rgt{:}),2);
end

if useCVmethod
    % Convert from degrees to radians
    angle_rad = unique_stimAngle / 180 * pi;
    
    % Take a vector average of all the orientations, weighted by average
    % response to each stimulus.
    osi = bsxfun(@times, average_responses, exp(1i*angle_rad));
    osi = sum(osi, 2);
    osi = bsxfun(@rdivide, osi, sum(average_responses,2));
    
    % Move from (-pi,pi] to [0,2pi) representation
    prefOri = mod(angle(osi), 2*pi);
    
    % Squash back down to [0,pi) if we were doing OSI instead of DSI
    if FLAG_OSI
        prefOri = prefOri / 2;
    end
    
    % Convert to degrees
    prefOri = radtodeg(prefOri);
    
    % Take the magnitude of the OSI as well
    prefSI = abs(osi);
else
    [pref_mag, pref_idx] = max(average_responses, [], 2);
    if FLAG_OSI
        unique_stimAngle = unique_stimAngle./2 ;% Squash back down to [0,135] if we were doing OSI instead of DSI
    end
    prefOri = unique_stimAngle(pref_idx);
    lowPref = pref_idx<=size(unique_stimAngle,2)/2;
    highPref = pref_idx>size(unique_stimAngle,2)/2;
    orth_idx = pref_idx;
    orth_idx = orth_idx + ((size(unique_stimAngle,2)/2)*lowPref);
    orth_idx = orth_idx - ((size(unique_stimAngle,2)/2)*highPref);  
    orth_mag = nan(size(pref_mag));
    for iDay = 1:size(average_responses,3)
        for iROI = 1:size(average_responses,1)
            for iAct = 1:size(average_responses,4)
                orth_mag(iROI,:,iDay,iAct) = average_responses(iROI,orth_idx(iROI,:,iDay,iAct),iDay,iAct);
            end
        end
    end
    prefSI = (pref_mag-orth_mag)./(pref_mag+orth_mag);
end

%retain nans in prefOri, otherwise prefOri will appear to be 0
prefOri(isnan(prefSI))=nan;

end