function ss = stimTimes2stimmap(stimTimes, tt)
% Label every point in time with the stimulus on screen

% Input handling ----------------------------------------------------------
assert(~iscell(stimTimes), 'stimTimes should be an array.');
assert(~iscell(tt)       , 'tt should be an array.');

% Main --------------------------------------------------------------------

ss = nan(size(tt));

nStim = numel(stimTimes);

% Handle case where there are no stimuli
if nStim==0
    ss(:) = 0;
    return;
end

% Set frames before the first stimulus to have no stim
ss(tt < stimTimes(1)) = 0;

% Handle labelling of each stimulus
for iStim=1:nStim
    li = tt >= stimTimes(iStim);
    if iStim<nStim
        li = li & tt < stimTimes(iStim+1);
    end
    ss(li) = iStim;
end

end
