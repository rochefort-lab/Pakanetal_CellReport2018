% Finds the indices at which locomotion begins and ends.
%
% Periods of locomotion are subdivided into 5 subperiods: before locomotion
% (A), immediately after locomotion onset (B), during locomotion (C),
% immediately before locomotion offset (D), and after locomotion offset
% (E).
% Subperiods A,B,D,E have a duration of 1 second each, whilst C has a
% duration as long as possible to extend from B to D.
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
% onsets : array [5, numLocomotionPeriods]
%     Indices of onsets of locomotion subperiods. The subperiods are
%     arranged along dimension 1, with rows 1-5 corresponding to A-E.
% offsets : array [5, numLocomotionPeriods]
%     Indices of offsets of locomotion subperiods. Every element in offsets
%     is the offset of the subperiod whose onset is in the corresponding
%     locotion in the onsets output.
%
% Notes
% -----
% To get the actual locomotion interval, use onsets(2,i):offsets(4,i),
% which runs from the beginning of B to the end of D.
%
% Illustration
% ------------
% Along dimension 1, of the returned onsets and offsets,
%     1: (A) before locomotion
%     2: (B) immediately after locomotion onset
%     3: (C) during locomotion
%     4: (D) immediately before locomotion offset
%     5: (E) after locomotion offset
%
%      | A | B |    C     | D | E |         vel
%      |   |   |          |   |   |          ^
%      |   |  _|__/\   /\_|   |   |          |
%      |   | / |    \_/   \__ |   |          |
%  _/\_|___|/  |          |  \|___|_____     +----> time
%  0220000033333333333333333333n00000000 - actlab(t)

%  0000000011111111111111111111000000000 :isloco
%  0000000+0000000000000000000-00000000  :diff(isloco)

function [onsets, offsets] = motion_onoffset(actlab, tt)

% -------------------------------------------------------------------------
% Parameters
duration         = 1.0 ; % seconds
max_dur_mergable = 0.5 ; % sec
min_still_prop   = 0.75; % should have at least 75% still to be valid A/E buffer

% -------------------------------------------------------------------------
% Sampling frequency is inverse of the difference in sampling times
%fs = 1./mean(diff(tt));
fs = 1./(tt(2)-tt(1));
duration_idx = round(duration * fs); % duration in number of sampling points

% Need to pad the action labels to find edges at the start/end of recording
pad_al = [NaN, actlab, NaN];

% We want to find locomotion ...
isloco     = pad_al==3;
% ... but it could have unmatched periods interspersed with it
islocoish  = pad_al==3 | pad_al==-1;
% and we want to find stationary ...
isstill    = pad_al==0;
% ... and it could have unmatched or twitch periods interspersed with it
isstillish = pad_al==0 | pad_al==-1 | pad_al==1 | pad_al==2;

% Find start and end of locomotion periods
loco_onset  = 1 + find(diff(isloco)== 1); % we add 1 so relative to pad_al start
loco_offset = 1 + find(diff(isloco)==-1); % we add 1 so relative to pad_al start

% Check whether we can merge any locomotion periods
ismergable = false(length(loco_onset)-1,1);
for iLoco = 1:length(loco_onset)-1
    % To be merged, we need to have a gap filled with "locoish" -1/1/3
    % and be short enough
    ismergable(iLoco) = ...
        all(islocoish(loco_offset(iLoco):loco_onset(iLoco+1))) & ...
        (loco_onset(iLoco+1)-loco_offset(iLoco)) / fs < max_dur_mergable;
end

% Merge all that we can
while any(ismergable)
    % Get the first in the list to merge
    iLoco = find(ismergable,1);
    % NB: The i-th element in ismergable is whether the i-th loco period
    % can be merged with its successor
    
    % Remove the end of the first in the chain
    loco_offset = loco_offset([1:(iLoco-1) (iLoco+1):end]);
    % Remove the start of the second in the chain
    loco_onset  = loco_onset([1:iLoco (iLoco+2):end]);
    % Remove the mergability status between the pair
    ismergable = ismergable([1:(iLoco-1) (iLoco+1):end]);
end

% Now check that the locomotion is long enough
% It has to be at least twice duration_idx, so that we can fit in
% a B and a D section
isLocoLong = (loco_offset - loco_onset) >= duration_idx * 2;

nLoco = length(loco_onset);

% init
onsets  = nan(5,nLoco);
offsets = nan(5,nLoco);

if nLoco==0
    % Bail now if there are no locomotion periods
    return;
end

% Approximately correct...
onsets(1,:)  = loco_onset - duration_idx    ; % A
onsets(2,:)  = loco_onset                   ; % B
onsets(3,:)  = loco_onset + duration_idx    ; % C
onsets(4,:)  = loco_offset- duration_idx    ; % D
onsets(5,:)  = loco_offset                  ; % E

offsets(1,:) = loco_onset                - 1; % A
offsets(2,:) = loco_onset + duration_idx - 1; % B
offsets(3,:) = loco_offset- duration_idx - 1; % C
offsets(4,:) = loco_offset               - 1; % D
offsets(5,:) = loco_offset+ duration_idx - 1; % E

% Remove datapoints which are out of range
% Check for out-of-bounds at the start
for iLoco=1:nLoco
    if onsets(1,iLoco) < 2
        % First A starts before data starts
        % Check whether B locomotion is from the beginning of the data
        if onsets(2,iLoco)<=2
            % Replace with C from the beginning
            onsets(3,iLoco) = 2;
            % Delete the on/offsets for B
            onsets(2,iLoco) = NaN;
            offsets(2,iLoco) = NaN;
            % Allow the loco is it is long enough for D
            % Doesn't need to be long enough for B
            isLocoLong(iLoco) = (loco_offset(iLoco) - loco_onset(iLoco)) >= duration_idx;
        end
        % Delete the on/offsets because the duration of A is insufficient
        onsets(1,iLoco) = NaN;
        offsets(1,iLoco) = NaN;
    end
    % Check for out-of-bounds at the end
    if offsets(5,iLoco) > length(pad_al)-1
        % Last E finishes after data ends
        % Check whether D locomotion is all the way to the end of the data
        if offsets(4,iLoco)>=length(pad_al)-1
            % Replace with C to the end
            offsets(3,iLoco) = length(pad_al)-1;
            % Delete the on/offsets for D
            onsets(4,iLoco) = NaN;
            offsets(4,iLoco) = NaN;
            % Allow the loco is it is long enough for D
            % Doesn't need to be long enough for B
            isLocoLong(iLoco) = (loco_offset(iLoco) - loco_onset(iLoco)) >= duration_idx;
        end
        % Delete the on/offsets because the duration of E is insufficient
        onsets(5,iLoco) = NaN;
        offsets(5,iLoco) = NaN;
    end
end

% Decided not to do this
% % % Shift back any A which ends with actlab==-1 or 1
% % for iLoco=1:nLoco
% %     % ...
% % end
% % % Shift foward any E which begins with actlab==-1 or 1
% % for iLoco=1:nLoco
% %     % ...
% % end
% % 
% % % Remove datapoints which are out of range
% % [onsets, offsets] = sanitize(onsets, offsets, length(pad_al));

% init
isPreValid  = false(1,nLoco);
isPostValid = false(1,nLoco);
for iLoco=1:nLoco
    if isnan(onsets(1,iLoco))
        % It's not valid because its out of range
        isPreValid(iLoco)  = false;
    else
        idx = onsets(1,iLoco):offsets(1,iLoco);
        % Check for sufficient A duration
        % With at least half actually time points actually stationary
%         isPreValid(iLoco)  = all(isstillish(idx)) & ...
%             sum(isstill(idx)) / length(idx) > min_still_prop
        isPreValid(iLoco)  = all(isstillish(idx));
    end
    if isnan(onsets(5,iLoco))
        % It's not valid because its out of range
        isPostValid(iLoco) = false;
    else
        idx = onsets(5,iLoco):offsets(5,iLoco);
        % Check for sufficient E duration
%         % With at least half actually time points actually stationary
%         isPostValid(iLoco) = all(isstillish(idx)) & ...
%             sum(isstill(idx)) / length(idx) > min_still_prop;
        isPostValid(iLoco)  = all(isstillish(idx));
    end
end

% Remove A and E which are too short
onsets( 1,~isPreValid ) = NaN;
offsets(1,~isPreValid ) = NaN;
onsets( 5,~isPostValid) = NaN;
offsets(5,~isPostValid) = NaN;

% Take off the padding
onsets  = onsets-1;
offsets = offsets-1;

% Drop ones which are not sufficiently long enough
onsets  = onsets(:,isLocoLong);
offsets = offsets(:,isLocoLong);

end
