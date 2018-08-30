function c = binsettings2correlation(binDur, binIntv)

fname = get_filename();

if ~exist(fname, 'file')
    generate_file();
end

% Get the look-up table
Cdat = load(fname);

% Interpolate our values against the look-up table
c = interp2(Cdat.binIntvs, Cdat.binDurs, Cdat.correlation, binIntv, binDur);

end


function fname = get_filename()

filename = 'binsample_correlations.mat';

% Find the correlation between consecutive samples

% Get the path of this function
% This is the path to getRepoDir.m
myFunPath = mfilename('fullpath');
% Get the folder part of the path
myFunPath = fileparts(myFunPath);

fname = fullfile(myFunPath, filename);

end


function generate_file()

disp('Creating bin sample correlation cached file');

Dat = load(fullfile(getRepoDir(), 'cache', ...
    'layer-sample_corr_nested-deltaf40Hz,ALL.mat'));

% bin duration, in seconds
Cdat.binDurs = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 1.9];
% bin interval, in seconds
Cdat.binIntvs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
% Average over all ROIs
Cdat.correlation = permute(nanmedian(Dat.data(:,:,:,1), 1), [2, 3, 1]);

disp(get_filename())

% Save file
save(get_filename(), '-struct', 'Cdat');

end