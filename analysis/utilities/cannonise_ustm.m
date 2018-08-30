function ustm_canon = cannonise_ustm(ustm)

% Remove the instance count from unique stimulus IDs
% We will let every stimulus be 
instanceMarkerStarts = regexp(ustm, ':\d*$', 'start');
assert(~any(cellfun(@isempty, instanceMarkerStarts)), ...
    'Some provided unique stim IDs do not have instance counts.');
ustm_canon = cell(size(ustm));
for iStm=1:numel(ustm_canon)
    ustm_canon{iStm} = ustm{iStm}(1:instanceMarkerStarts{iStm}-1);
end

end