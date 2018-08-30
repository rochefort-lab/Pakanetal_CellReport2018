% Downsample the data set for a single recording.
%
% The inputs are always decimated along dimension 2. The inputs can must be
% 2D or 3D, and may have any size along dimensions 1 and 3.
%
% Inputs
% ------
% target_fs : double
%     The sampling frequency at which to resample and return the data.
% signal : array [*, numTimePoints, ...]
%     Raw signal to be decimated, with the signal lowpass filtered and
%     downsampled.
% tt : array [*, numTimePoints, ...]
%     Times of sampling points (in seconds). The resampling factor (R) is
%     determined automatically from the original sampling frequency
%     implicitly given by tt and the supplied target_fs frequency.
% actlab : array [*, numTimePoints, ...] (optional)
%     Action labels, to be downsampled without lowpass filtering.
% *varargin : arrays each [*, numTimePoints, ...] (optional)
%     Additional inputs to be decimated, all with lowpass filtering before
%     dowsampling.
%
% Outputs
% -------
% signal : array [*, numTimePoints/R, ...]
%     Decimated signal.
% tt : array [*, numTimePoints/R, ...]
%     Times of sampling points.
% actlab : array [*, numTimePoints/R, ...] (optional)
%     Action labels, downsampled without lowpass filtering.
% *varargout : arrays each [*, numTimePoints/R, ...] (optional)
%     Additional outputs also decimated from varargin.
%
% See also decimate2.

function [signal, tt, varargout] = decimate_rec(target_fs, signal, tt, varargin)

% Check inputs
if ndims(signal)>3
    error('Signal should be two- or three-dimensional');
end
if ndims(tt)>3
    error('Timepoints input should be two- or three-dimensional');
end

% Measure input frequency
dt = diff(tt, 1, 2);
fs = 1./mean(dt(:));

% Check the resampling interval
R = fs./target_fs;

if abs(R - round(R)) > 1e-10
    warning('Resampling rate is not integer')
end
R = round(R);

if R==1
    % No need to resample
    varargout = varargin;
    return;
elseif R<1
    error('Can''t sample at %.2fHz - input only at %.2fHz', target_fs, fs);
end

% Decimate every time-signal
signal = decimate2(signal, R, 2, 'FIR');
if ischar(varargin{1,4}) || size(varargin{1,4},2) == 1
    for iArg=2:numel(varargin)-1
    % ped, vel, acc
    varargin{iArg} = decimate2(varargin{iArg}, R, 2, 'FIR');
    end
else
    for iArg=2:numel(varargin)
    % ped, vel, acc
    varargin{iArg} = decimate2(varargin{iArg}, R, 2, 'FIR');
    end

end
% Downsample time and actionlabels
tt = tt(:, 1:R:end, :);
if ~isempty(varargin)
    % actlab
    varargin{1} = varargin{1}(:,1:R:end,:);
end

% Define output
varargout = varargin;

end
