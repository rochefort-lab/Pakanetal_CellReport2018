% Compute LMI from samples for each behavioural condition.
%
% Inputs
% ------
% signal : array [numRois, numSamples]
%     The signal of interest.
% actlab : array [1, numSamples]
%     Which action is being performed during each sample, with stationary
%     samples denoted with 0 and locomotion with 3.
% ss_string : cell array [1, numSamples]
%     Which stimulus is being presented during each sample, in cannonical
%     form. The first character of each stimulus string must be the
%     category of the stimulus.
% discard_uniform : bool, optional
%     If true, discard any samples from uniform gratings. Default is
%     `false`.
%
% Outputs
% -------
% lmi : array [numRois, 1]
%     Unbiased estimate of locomotion modulation index across all stimuli.
%     Each stimulus category (such as static, forward, uniform) is weighted
%     evenly. Within each category, each stimulus is weighted evenly to get
%     an unbiased estimate for that category.
%
% See also: binnedrec2samples.

function lmi_avg = samples2lmi(signal, actlab, stm, discard_uniform)

% Input handling
if nargin<4
    discard_uniform = false;
end

% Discard uniform stimuli, if requested
if discard_uniform
    is_uniform = strncmp(stm, 'U', 1);
    signal = signal(:, ~is_uniform);
    actlab = actlab(:, ~is_uniform);
    stm = stm(:, ~is_uniform);
end

[uu, ~, uu_idx] = unique(stm);  % list every unique stimulus
uu_idx = uu_idx';  % Fix vector orientation
aa = [0, 3];  % only care about stationary or locomotion samples
mm = nan(size(signal,1), numel(uu), numel(aa));  % for mean of signal
nn = zeros(1, numel(uu), numel(aa));  % count of {action, stim} occurances

% Take the mean of signal for every behavioural condition (action and stim)
for iStm=1:numel(uu)
    is_stm = uu_idx==iStm;
    for iAct=1:numel(aa)
        is_act = actlab==aa(iAct);
        is_ours = is_stm & is_act;
        mm(:, iStm, iAct) = nanmean(signal(:, is_ours), 2);
        nn(:, iStm, iAct) = sum(is_ours);
    end
end

% if any values are negative, set them to nan
mm(mm<0)=nan;

% LMI is the relative difference between locomotion and stationary
lmi = (mm(:, :, 2) - mm(:, :, 1)) ./ (mm(:, :, 2) + mm(:, :, 1));
% Set any datapoints which had no samples for one or the other to be NaN
lmi(:, any(nn==0, 3)) = NaN;

% Assume the first letter of the ustm is the stimulus category
ustm_first_letters = cellfun(@(x) x(1), uu);
% We don't care if stimulus is moving forward or backward, just that its
% moving, so replace every 'B' with 'F'
ustm_first_letters(ustm_first_letters=='B') = 'F';
unique_first_letters = unique(ustm_first_letters);

% Take the average LMI for each stimulus category
lmi_by_cat = nan(size(lmi,1), numel(unique_first_letters));
for iCat=1:numel(unique_first_letters)
    is_cat = ustm_first_letters == unique_first_letters(iCat);
    lmi_by_cat(:, iCat) = nanmean(lmi(:, is_cat), 2);
end

% Take the average over all categories
lmi_avg = nanmean(lmi_by_cat, 2);

end