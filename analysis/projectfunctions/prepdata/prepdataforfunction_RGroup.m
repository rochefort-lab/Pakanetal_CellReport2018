function [varargout] = prepdataforfunction_RGroup(nameOfFunction, varargin)

List_func_do_VRdata2tracks = {...
    'binnedVRData','gainModVRData','pairedTestGratVSrwdZone',...
    'tempMatchVR_RZ','tempMatchVR_HitMiss','slowMedfastTrials','mancova'};


if ismember(nameOfFunction, List_func_do_VRdata2tracks)
    % *********************************************************************
    % -- GROUP OF FUNCTIONS CALLING VRdata2tracks
    % *********************************************************************
    
    signal_df = varargin{1};
    time_samples = varargin{2};
    VRlab = varargin{3};
    disp = varargin{4};
    vel = varargin{5};
    lagTime = varargin{6};
    idxArgin_suppl = 8;
    % DEFAULT PARAMS
    dataTrialsOnly = []; % if true, only returns trials with data (remove nan trials of other tracks or categories)
    useHitsOnly = false; % if true, only successful trials are used (licked in reward zone)
    useMissOnly = false;
    
    switch nameOfFunction
        
        case 'binnedVRData'
            dataTrialsOnly = 1; % 1 = get data trials only
            numBin_max = 96; % Max number of bins
            
            if nargin >= idxArgin_suppl && ~isempty(varargin{idxArgin_suppl})
                switch varargin{idxArgin_suppl}
                    case 'useHitsOnly'
                        useHitsOnly = true;
                    case 'useMissOnly'
                        useMissOnly = true;
                end
            end
            
        case {'gainModVRData' , 'mancova'}
            dataTrialsOnly = true; % 1 = get data trials only
            
        case 'tempMatchVR_RZ'
            useHitsOnly = true;
            
    end
    
    % Organize by track/distance etc
    [~, dFTrials, ~, meandF, ~, ~, hitTrial, velTrials, trialTimes] = ...
        VRdata2tracks(signal_df, time_samples, VRlab, disp, vel, lagTime, [], [], [], dataTrialsOnly, useHitsOnly, useMissOnly);
    % Note:
    % dFTrials:  cell array [numTracks, 1]
    %   each cell is: array [nROIs, numBins, numTrials]
    %       Binned fluorescence dF/F signal for each trial of the track
    
    % Check empty tracks
    nonemptyTracks = ~ cellfun(@isempty, dFTrials);
    
    % Get number of tracks
    nTracks = length(dFTrials);
    
    % Set output depending on function requested
    switch nameOfFunction
        
        case 'binnedVRData'
            varargout{1} = meandF;
            varargout{2} = numBin_max;
            
        case 'gainModVRData'
            takeTrials = 5; % take every the nth trial (starting with trial trialNumStart)
            trialNumStart = 4; % take every the nth (defined by takeTrials) trial (i.e starting with trial 4)
            numBin_max = 96; % Max number of bins
            
            % === Separate the gain trials and others
            dataGainTrials = dFTrials; % Initialise
            dataNonGainTrials = cell(nTracks, 1);
            for iTrack = 1:2:nTracks
                if nonemptyTracks(iTrack)
                    numTrials = size(dFTrials{iTrack},3);
                    Idx_gainTrials = false(1,numTrials);
                    Idx_gainTrials(trialNumStart:takeTrials:numTrials) = true;
                    dataGainTrials{iTrack} = dFTrials{iTrack}(:, :, Idx_gainTrials);
                    dataNonGainTrials{iTrack} = dFTrials{iTrack}(:, :, ~Idx_gainTrials);
                end
            end
            
            varargout{1} = dataGainTrials;
            varargout{2} = dataNonGainTrials;
            varargout{3} = numBin_max;
            
        case 'pairedTestGratVSrwdZone'
            % Note: RZ starts 16 bins (40cm) from end; 120cm length track RZ starts at bin 32
            % To have same amt of data (10 bins)
            dataRangepre = [30,20]; % from end of trial; Gratings: 10 bins total (25cm), from end [30,20] = -35 to -10 cm from reward zone onset.
            dataRangepost = [16,6]; % from end of trial; Reward Zone: 10 bins total (25cm), from end [16,6] = 0 to +25 cm from reward zone onset.
            
            RZdatapre = cell(nTracks, 1);
            RZdatapost = cell(nTracks, 1);
            
            RZdatapre(nonemptyTracks) = cellfun(@squeeze, (cellfun(@(x) nanmean(x(:,end-dataRangepre(1):end-dataRangepre(2),:),2), dFTrials(nonemptyTracks), 'UniformOutput', false)),'UniformOutput', false);
            RZdatapost(nonemptyTracks) = cellfun(@squeeze, (cellfun(@(x) nanmean(x(:,end-dataRangepost(1):end-dataRangepost(2),:),2), dFTrials(nonemptyTracks), 'UniformOutput', false)),'UniformOutput', false);
            
            varargout{1} = RZdatapre;
            varargout{2} = RZdatapost;
                
        case 'tempMatchVR_HitMiss'
            
            dataDf = cell(nTracks, 1);
            dataDf(nonemptyTracks) = cellfun(@squeeze, (cellfun(@(x) nanmean(x,2), dFTrials(nonemptyTracks), 'UniformOutput', false)),'UniformOutput', false);
            
            % Separate the sucessful (Hit) and miss trials
            dataDfHit = cell(nTracks, 1);
            dataDfMiss = cell(nTracks, 1);
            for iTrack = 1:nTracks
                if nonemptyTracks(iTrack)
                    dataDfHit{iTrack} = dataDf{iTrack}(:, hitTrial{iTrack});
                    dataDfMiss{iTrack} = dataDf{iTrack}(:, ~hitTrial{iTrack});
                end
            end
            
            varargout{1} = dataDfHit;
            varargout{2} = dataDfMiss;
            
        case 'tempMatchVR_RZ_CuedUncued'
            dataRange = [18,1];  % from end of trial, RZ starts 16 bins from end; data taken from 5cm before and throughout the reward zone.
            trackType_Idx = {[1,2]}; % compare visually cued (track 1) to uncued (track 2)
            
            RZdataDf = cell(nTracks, 1);
            
            RZdataDf(nonemptyTracks) = cellfun(@squeeze, (cellfun(@(x) nanmean(x(:,end-dataRange(1):end-dataRange(2),:),2), dFTrials(nonemptyTracks), 'UniformOutput', false)),'UniformOutput', false);
            
            varargout{1} = RZdataDf;
            varargout{2} = trackType_Idx;
            
            
            
        case 'slowMedfastTrials'
            dataRange = [20,12]; % from end of trial
            
            dataDf_RZone_slow = cell(nTracks, 1);
            dataDf_RZone_med = cell(nTracks, 1);
            dataDf_RZone_fast = cell(nTracks, 1);
            % cell array [numTracks, 1]
            %   each cell is: array [nROIs, numTrials]
            
            % Get data surrounding reward zone (bin 16 from end, 8 bins surrounding  = 10 cm total)
            for iTrack = 1:nTracks
                is_slow_trial = trialTimes{iTrack} >= prctile(trialTimes{iTrack}, 75); % slowest 25%   size: [numTrials, 1]
                dFTrial_slow = dFTrials{iTrack}(:,:,is_slow_trial); % dF{iTrack} is [nROIs, numBins, numTrials]
                dataDf_RZone_slow{iTrack} = squeeze( nanmean(dFTrial_slow(:,end-dataRange(1):end-dataRange(2),:),2) );
                
                is_med_trial = (trialTimes{iTrack} >= prctile(trialTimes{iTrack}, 25)) & (trialTimes{iTrack} <= prctile(trialTimes{iTrack}, 75));
                dFTrial_med = dFTrials{iTrack}(:,:,is_med_trial);% cellfun(@(x, y) y(:,:,logical(x(1,1,:))), fastTrials, dFTrials, 'UniformOutput', false);
                dataDf_RZone_med{iTrack} = squeeze( nanmean(dFTrial_med(:,end-dataRange(1):end-dataRange(2),:),2) );                
                
                is_fast_trial = trialTimes{iTrack} <= prctile(trialTimes{iTrack}, 25); % fastest 25% size: [numTrials, 1]
                dFTrial_fast = dFTrials{iTrack}(:,:,is_fast_trial);
                dataDf_RZone_fast{iTrack} = squeeze( nanmean(dFTrial_fast(:,end-dataRange(1):end-dataRange(2),:),2) );
            end        
            
            varargout{1} = trialTimes;
            varargout{2} = dataDf_RZone_slow;
            varargout{3} = dataDf_RZone_med;
            varargout{4} = dataDf_RZone_fast;
            
        case 'mancova'
            is_shuffle_case = nargin == 8 && ~isempty(varargin{idxArgin_suppl}) && varargin{idxArgin_suppl};
            
            dataRange = [20,12]; % from end of trial

            varargout{1} = trialTimes;
            varargout{2} = hitTrial;
            
            if is_shuffle_case
                varargout{3} = cellfun(@squeeze, (cellfun(@(x) permute(x, [1,3,2]), dFTrials, 'UniformOutput', false)),'UniformOutput', false);
            else
                dataDf_RZ = cellfun(@squeeze, (cellfun(@(x) nanmean(x(:,end-dataRange(1):end-dataRange(2),:),2), dFTrials, 'UniformOutput', false)),'UniformOutput', false);
                varargout{3} = dataDf_RZ;
            end
            
            dataVel_RZ = cellfun(@squeeze, (cellfun(@(x) nanmean(x(:,end-dataRange(1):end-dataRange(2),:),2), velTrials, 'UniformOutput', false)),'UniformOutput', false);
            varargout{4} = dataVel_RZ;
    end
    
else
    % *********************************************************************
    % -- VARIOUS FUNCTIONS USING VR data
    % *********************************************************************
    
    
    
    switch nameOfFunction
        
        case 'rwdOnsetData'
            
            signal_df = varargin{1};
            time_samples = varargin{2};
            VRlab = varargin{3}; % note: all fields are [1, numSamples]
            
            if length(varargin)>3 && ~isempty(varargin{4})
                track_distance = varargin{4};
                doRewardZoneOnset = true;
            else
                doRewardZoneOnset = false;
            end
            
            % DEFAULT PARAMS of tracks
            [numTracks, TracksParameters] = get_tracks_info();
            
            avgHz = 1/(mean(diff(time_samples)));
            window_samples = floor(avgHz*2); % (in frames): avgHz = 1000ms window; get 2sec before and after Rwd Onset
            window_time = time_samples(1:window_samples*2);
            nROIs = size(signal_df,1);
            
            dataDf_RewardWindow = cell(1,2);
            numSamples = length(VRlab.trackNum);
            list_indices = 1:numSamples;
            
            % If detect reward zone rather than reward onset, find for each
            % trial the index of the onset of the reward zone depending on
            % the type of reward (correct or default)
            % NOTE: skip last trial in case it is incomplete
            if doRewardZoneOnset
                numTrials = max(VRlab.trialNum)-1;
                distOnsetEarly = nan(1,numTrials);
                distOnsetLate = nan(1,numTrials);
                for iTrial = 1:numTrials
                    is_sample_of_trial = VRlab.trialNum == iTrial;
                    idx_sample_startTrial = find(is_sample_of_trial,1);
                    track_num = VRlab.trackNum(idx_sample_startTrial);
                    
                    idx_sample_onset = idx_sample_startTrial + ...
                        find(track_distance(is_sample_of_trial)>TracksParameters(track_num).rewardZone,1);
                    % Check that event was not too soon or too late
                    is_possible_Index = ( idx_sample_onset >= (window_samples*2 +1) ) & ( (idx_sample_onset+window_samples) <= numSamples );
                    
                    if is_possible_Index
                        if any(VRlab.reward(is_sample_of_trial)==1) % Find if reward hit
                            distOnsetEarly(iTrial) = idx_sample_onset;
                            
                        else % Default reward
                            distOnsetLate(iTrial) = idx_sample_onset;
                        end
                    end
                end
            else
                is_possible_Index = ( list_indices >= (window_samples*2 +1) ) & ( (list_indices+window_samples) <= numSamples );
            end
            
            
            for iRwdType = 1:2
                % Find reward events
                if doRewardZoneOnset
                    if iRwdType==1
                        list_indices_selected = distOnsetEarly(~isnan(distOnsetEarly));
                    else
                        list_indices_selected = distOnsetLate(~isnan(distOnsetLate));
                    end
                    if isempty(list_indices_selected)
                        continue;
                    end
                else
                    is_reward_event = VRlab.reward==iRwdType;
                    
                    if ~any(is_reward_event)
                        continue
                    end
                    
                    % Make sure that an event was not too soon or too late
                    is_reward_event = is_reward_event & is_possible_Index ;
                    list_indices_selected = list_indices(is_reward_event);
                end
                
                numSamples_events = length(list_indices_selected);
                
                eventAvg = nan(nROIs, window_samples*2, numSamples_events);
                for iEvent = 1:numSamples_events
                    selected_samples = list_indices_selected(iEvent)-(window_samples-1) : list_indices_selected(iEvent)+window_samples ;
                    % Take mean data across trials from window after events (onset | dataWindow)
                    eventAvg(:,:,iEvent) = signal_df(:,selected_samples);
                end
                
                dataDf_RewardWindow{iRwdType} = cell(1,numTracks);
                for iTrack = 1:numTracks
                    is_event_of_track = VRlab.trackNum(list_indices_selected) == iTrack;
                    dataDf_RewardWindow{iRwdType}{iTrack} = eventAvg(:,:,is_event_of_track);
                end
            end
            
            varargout{1} = dataDf_RewardWindow;
            varargout{2} = window_samples;
            varargout{3} = window_time;
            varargout{4} = numTracks;
            
            
        case 'successRateSMI_VR'
            
            VRlab = varargin{1}; % note: all fields are [1, numSamples]
            track_distance = varargin{2};
            
            % DEFAULT PARAMS of tracks
            [numTracks, TracksParameters] = get_tracks_info();
            
            success_rate = nan(1, numTracks);
            dataReward = cell(1,numTracks);
            dataLicks = cell(1,numTracks);
            dataDistance = cell(1,numTracks);
            
            % Calculate success rate and get reward & lick data
            for iTrack = 1:numTracks
                track_samples = VRlab.trackNum == iTrack;
                if ~any(track_samples)
                    continue
                end
                
                list_trialNumbers_track = unique(VRlab.trialNum(track_samples));
                numTrials = length(list_trialNumbers_track);
                dataReward{iTrack} = cell(1, numTrials);
                dataLicks{iTrack} = cell(1, numTrials);
                dataDistance{iTrack} = cell(1, numTrials);
                numhitTrials = 0;
                for iTrial = 1:numTrials
                    selected_trial = list_trialNumbers_track(iTrial);
                    trial_track_samples = track_samples & (VRlab.trialNum == selected_trial);
                    if any(VRlab.reward(trial_track_samples)==1)
                        numhitTrials = numhitTrials+1;
                    end
                    
                    dataReward{iTrack}{iTrial} = VRlab.reward(trial_track_samples);
                    dataLicks{iTrack}{iTrial} = VRlab.licks(trial_track_samples);
                    dataDistance{iTrack}{iTrial} = track_distance(trial_track_samples);
                end
                
                success_rate(iTrack) = numhitTrials/numTrials; % percent success rate, num hits, num misses
            end
            
            
            
            % Output
            varargout{1} = success_rate;
            varargout{2} = dataReward;
            varargout{3} = dataLicks;
            varargout{4} = dataDistance;
            varargout{5} = TracksParameters;
            
            
            
    end
    
end



end

%__________________________________________________________________________
% ___ FUNCTIONS ___________________________________________________________


function [numTracks, TracksParameters] = get_tracks_info()


numTracks = 8;
% Tracks parameters
TracksParameters = struct('rewardZone', cell(1,numTracks), 'rewardDefault', cell(1,numTracks));
[TracksParameters([1,2]).rewardZone] = deal(80);
[TracksParameters([1,2]).rewardDefault] = deal(100);
[TracksParameters([3,5]).rewardZone] = deal(120);
[TracksParameters([3,5]).rewardDefault] = deal(140);
[TracksParameters([7,8]).rewardZone] = deal(200);
[TracksParameters([7,8]).rewardDefault] = deal(220);

end
