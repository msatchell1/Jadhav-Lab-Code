% Script to calculate the ripple occurance start and end time structures
% using the new methods I have written. Results are saved to the animal's
% folder.
% Michael Satchell 08/28/2023


%% Load data

clearvars;

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
load_rats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','pos01'
filetypes = {'spikes01','tetinfo','pos01'};

C_alldata = {}; % Cell array to hold data for all rats. If multiple filetypes 
% are loaded, each row holds a different file type, ordered in the same
% order as the elements in filetypes.

disp("Loading new animal data... ")
for r = 1:length(load_rats)
    fprintf("Loading animal: %s \n",load_rats{r})

    for ft = 1:length(filetypes)    

        short_name = load_rats{r};
        chop_idx = strfind(load_rats{r},'_') - 1;
        if ~isempty(chop_idx)
            short_name = load_rats{r}(1:chop_idx); % Gets the first characters of the rat's name before an '_'.
            % So far this is only needed for ER1_NEW to remove the '_NEW'.
        end

        % Does not load from EEG folder
        File_dir = dir(data_dir+"/"+load_rats(r)+'_direct'+"/"+short_name+filetypes{ft}+"*");
    
        if isempty(File_dir)
            error("%s file does not exist for animal: %s \n",filetypes{ft},load_rats{r})
        elseif length(File_dir) > 1
            error("More than one file detected when searching for: %s, in animal: %s \n" + ...
                "Change names in filetypes to be more specific.", filetypes{ft},load_rats{r});
        else
            file = struct2cell(load(string(fullfile(File_dir.folder, File_dir.name)))); % load data
            file = file{:};
            C_alldata{ft,r} = file{1,1};
            fprintf("       Loaded file: %s  \n",  File_dir.name)
        end
    end
end

C_alldata = clip_17_epochs(C_alldata); % removes extra epoch data.

spikes_idx = find(contains(filetypes,'spikes01')); 
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end

pos_idx = find(contains(filetypes,'pos01'));
if isempty(pos_idx)
    error("pos01 data must be loaded to run this analysis.")
end

ti_idx = find(contains(filetypes,'tetinfo'));
if isempty(ti_idx)
    error("tetinfo data must be loaded to run this analysis.")
end

C_allspikes = C_alldata(spikes_idx,:);
C_runstate = create_runstate(C_alldata(pos_idx,:));
C_stillstate = create_stillstate(C_alldata(pos_idx,:));
C_alltetinfo = C_alldata(ti_idx,:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;


%% Create new rippletime files with improved SWR detection

% All riptets for each animal:

% For ZT2
% riptetsStr = {'10','11','12','13','14','15','16','17', '18','22','23', '19','20','21','24','25',...
%     '26','27','28','29','34','36', '31','32','33','35'};

% For ER1_NEW
% riptetsStr = {'26','25','24','23', '11','07','21','08','15','13', '14'};

% For KL8
% riptetsStr = {'09','10', '25','24','15','21','23','22'}; 

% For BG1
% riptetsStr = {'12','10','07','21','22','25'}; 

% For JS14
% riptetsStr = {'12','06','08','09','10','15','26'}; 

% For JS15
% riptetsStr = {'05','06','07','08','09','11','21','24','25'};

% For JS17
% riptetsStr = {'06','05','24','22','11','12'};

% For JS21
% riptetsStr = {'05','06','25','07', '08'};

% For JS34
% riptetsStr = {'06','07', '11','12', '22','24','25','09'};


% A cell array that holds riptets in their plotting order for all rats.
% NOTE: Assumed order of rats in 'load_rats' is: 
% {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'}
% IF load_rats ORDER CHANGES, C_riptets MUST BE ADJUSTED TO MATCH.
C_riptets = {
{'10','11','12','13','14','15','16','17', '18','22','23', '19','20','21','24','25','26','27','28','29','34','36', '31','32','33','35'},...
{'26','25','24','23', '11','07','21','08','15','13', '14'},...
{'09','10', '25','24','15','21','23','22'},...
{'12','10','07','21','22','25'},...
{'12','06','08','09','10','15','26'},...
{'05','06','07','08','09','11','21','24','25'},...
{'06','05','24','22','11','12'},...
{'05','06','25','07', '08'},...
{'06','07', '11','12', '22','24','25','09'},...
};

% % I want to save a struct with similar fields to that of Justin's
% % rippletimes struct. However, I want to have multiple versions of ripple
% % detection, so I will store multiple structs in a cell array for every
% % epoch in rippletimes. The first one will be Justin's rippletime data,
% % followed by my own ripple detection algorithm results. This will then be
% % saved to the animal's file, just like Justin's rippletimes. Then, when I load the cell
% % array I will have access to any of the different ripple detection
% % algorithms I applied. Each struct will have fields indicating the parameters
% % or methods used for detection (1/0 for each field might be best).
% % ripMethods will be the cell array.
% ripMethods = {};
% 
% % Ripple detection method names.
% methodNames = {"Justin"};

% C_allripmeth = cell(1,length(load_rats));
rStr = "ER1_NEW"; % rat name
r = find(contains(load_rats,rStr));
% for r = 1:length(load_rats)
    % rStr = load_rats{r};
    short_name = load_rats{r};
    chop_idx = strfind(load_rats{r},'_') - 1;
    if ~isempty(chop_idx)
        short_name = load_rats{r}(1:chop_idx); % Gets the first characters of the rat's name before an '_'.
        % So far this is only needed for ER1_NEW to remove the '_NEW'.
    end
    
    rippletimes = cell(1,17); % To hold rippletime structures for each epoch.
    for e = 1:17
    
        if e < 10 
            eStr = "0"+num2str(e); % epoch number
        else
            eStr = num2str(e);
        end
        
        
        for rt = 1:size(C_riptets{1,r},2)
        
            tetStr = C_riptets{1,r}{1,rt};
            tet = str2double(tetStr);
        
            eegref = load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/%s_direct/EEG/%seegref01-%s-%s.mat",rStr,short_name,eStr,tetStr));
            eegripple = load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/%s_direct/EEG/%sripple01-%s-%s.mat",rStr,short_name,eStr,tetStr));
            
            refData = eegref.eegref{1,1}{1,e}{1,tet};
            timeRange = (refData.starttime : 1/refData.samprate : refData.endtime)';
            refV = refData.data/refData.voltage_scaling; % Remove scaling factor
            refV = refV/1000; % Change from microvolts to millivolts.
            
            ripData = eegripple.ripple{1,1}{1,e}{1,tet};
            ripAmp = ripData.data(:,1)/ripData.voltage_scaling; % Ripple filtered LFP amplitude.
            ripAmp = ripAmp/1000; % Change from microvolts to millivolts.
            ripEnv = ripData.data(:,3); % Ripple filter envelope.
            
            % Zscore the envelope of the ripple filter to determine significant
            % ripples.
            zEnv = zscore(double(ripEnv));
        
            if rt == 1 % define matrix to hold envelope traces. Envelope data should be the same length for all tetrodes.
                riptetszEnvs = zeros([size(zEnv,1), size(C_riptets{1,r},2)]);
            end
        
            riptetszEnvs(:,rt) = zEnv; % add envelope data
        
        end
        
        gwinDur = 100; % Duration of gaussian window in ms.
        gwinlen = (refData.samprate/1000)*gwinDur + 1; % window length in time steps.
        gwin = gausswin(gwinlen);
        
        % Use neuron spike data only from CA1 neurons
        nrnTets = []; % Tetrodes to plot neurons from.
        % Sort out CA1 tetrodes from PFC ones.
        for tet = 1:size(C_alltetinfo{1,r}{1,e},2)
            tetinfo = C_alltetinfo{1,r}{1,e}{1,tet};
            if isfield(tetinfo,"area") && strcmp(tetinfo.area,"CA1")
                nrnTets(end+1) = tet;
            end
        end
        
        
        nrnCount = 0;
        nrnsGtrace = []; % Holds gaussian trace for all neurons.
        for i = 1:length(nrnTets)
            tet = nrnTets(i);
            if ~isempty(C_allspikes{1,r}{1,e}{1,tet})
                nrns = [C_allspikes{1,r}{1,e}{1,tet}(:)]; % cell array of structs, each row a neuron.
                for j = 1:length(nrns)
                    if isfield(nrns{j}, "meanrate") % This also checks that the struct exists
                        nrn = nrns{j};
                        meanRate = nrn.meanrate;
                        if meanRate < 7
                            spikeTimes = nrn.data(:,1);
            
                            isSpike = ismembertol(timeRange,spikeTimes, 1/(2*refData.samprate), DataScale=1);
                            gtrace = conv(isSpike,gwin,"same");
                            cutgtrace = gtrace > 1; % where to cut off summed peaks
                            gtrace(cutgtrace) = 1; % reassign these values to 1.
            
                            nrnCount = nrnCount + 1;
                           
                            nrnsGtrace(:,nrnCount) = gtrace;  % add nrn gaussian trace to matrix
                        end
                    end
                end
            end
        end


        
        % Sum the gaussian traces from all neurons and zscore it.
        sumGtrace = sum(nrnsGtrace,2);
        zGtrace = zscore(sumGtrace);
        
        % Calculate zscore of sum of envelope zscores of all tetrodes and plot.
        sumzEnv = sum(riptetszEnvs,2);
        popEnvzscore = zscore(sumzEnv); % Zscore of all tetrodes
        
        cumEnvG = zGtrace+popEnvzscore; % The cumulative zscore of both the spike data gaussians and LFP envelope.
        % zCum = zscore(zGtrace+popEnvzscore); % z-score the combined mulitunit and tetrode zscore traces.
        
        peakThr = 6; % Number of standard deviations above the mean to detect peak
        peakSep = 0;  % Minimum time (in ms) between peaks.
        peakWidth = 10; % Minimum width (in ms) that a peak must have.
        peakProm = 6; % Minimum height (in stds) above neighboring signal a peak must have.
        
        % Use the findpeaks function to find where the zscore passes above
        % threshold. Note that when the sample rate is provided, peak locations
        % and widths are returned in units of time.
        [pkHeights,pkTimes,pkWidths,pkProms] = findpeaks(cumEnvG,timeRange, MinPeakHeight=peakThr,...
            MinPeakDistance=peakSep/1000, MinPeakWidth=peakWidth/1000, MinPeakProminence=peakProm);
        
        % Peaks that occur during outside of still periods are removed. These are mostly
        % large amplitude noise events during behavior from the rat bumping
        % into things.
        stillStarttimes = C_stillstate{1,r}{1,e}.starttime;
        stillEndtimes = C_stillstate{1,r}{1,e}.endtime;
        indsToRemove = [];
        for p = 1:length(pkTimes)
            pkt = pkTimes(p);
            if all(~(pkt >= stillStarttimes & pkt <= stillEndtimes))
                % fprintf("Peak at time %d removed",pkt)
                indsToRemove(end+1) = p;
            end
        end  
        pkHeights(indsToRemove) = []; pkTimes(indsToRemove) = []; 
        pkWidths(indsToRemove) = []; pkProms(indsToRemove) = []; % Removes peaks
    
        % Now that the peaks have been found, it is time to determine the window
        % edges for rippletimes. To do this, the times at which cumEnvG crosses below 0 will
        % define the edges of the window for each peak. If multiple peaks occur
        % within the same window, only one rippletime window will be recorded. The
        % peak data should be saved separately.
        
        % Indices where cumulative zscore < 0.
        isltz = cumEnvG < 0;
        % If the zscore dips down very briefly below zero then comes back up I may
        % still want to include the whole thing as a single ripple, so a minimum
        % dip size is used to remove groups of indices in isltz below that
        % size.
        minDipDur = 5; % Minimum time (ms) required for zscore to be below 0 to exit ripple window.
        minDipLen = ceil((refData.samprate/1000)*minDipDur);
        
        % I need to group the indicies in ltzInds, but without being able to use
        % bwlabel (which I need image processing toolbox for, and can't install for
        % some reason), all the group sorting only sorts by value. So I need to
        % assign a different value to each group. For now I will do this with a for
        % loop...
        ltzLabel = -1;
        gtzLabel = 2;
        groupLabels = zeros(size(isltz));
        for i = 1:length(isltz)-1
            if isltz(i) == 0 % if zscore greater than zero
                groupLabels(i) = gtzLabel;
                if isltz(i+1) == 1
                    gtzLabel = gtzLabel + 1;
                end
            elseif isltz(i) == 1 % if zscore less than zero
                groupLabels(i) = ltzLabel;
                if isltz(i+1) == 0
                    ltzLabel = ltzLabel - 1;
                end
            end
        end
        
        % gcounts counts the number of time steps in each group and keeps the label
        % associated with those groups in gIDs.
        [gcounts, gIDs] = groupcounts(groupLabels);
        % Converts groups that are too small to be as if they are greater than zero
        % (i.e., assigned to value 0 in isltz).
        for g = 1:length(gcounts)
            if gcounts(g) < minDipLen
                isltz(groupLabels == gIDs(g)) = 0;
            end
        end
        
        starttimes = zeros(size(pkTimes));
        endtimes = zeros(size(pkTimes));
        % Now the edge determination can happen.
        for p = 1:length(pkTimes)
            pkt = pkTimes(p); % time in seconds of peak
            % pkidx = find(timeRange == pkt); % index in timeRange of peak
            isltpk_isltz = timeRange < pkt & isltz == 1; % times less than peak time and with zscore less than zero.
            isgtpk_isltz = timeRange > pkt & isltz == 1; % times greater than peak time and with zscore less than zero.
            starttimeidx = find(isltpk_isltz,1,"last");
            endtimeidx = find(isgtpk_isltz,1,"first");
            % Actually record start and end times around the peak.
            starttimes(p) = timeRange(starttimeidx);
            endtimes(p) = timeRange(endtimeidx);
        end
        
        % Remove duplicate times for ripples with multiple peaks.
        starttimes = unique(starttimes,"stable");
        endtimes = unique(endtimes,"stable");
        
        % make isRippleTime vector
        isRippleTime = zeros(size(timeRange));
        for t = 1:size(starttimes,1)
            isRippleTime(timeRange >= starttimes(t) & timeRange <= endtimes(t)) = 1; % For a single ripple
        end
        
        
        
        % S_rippletimes.methodName
        S_rippletimes.starttime = starttimes;
        S_rippletimes.endtime = endtimes;
        S_rippletimes.isRippleTime = isRippleTime;
        S_rippletimes.timeVec = timeRange; % in seconds
        S_rippletimes.total_duration = sum(endtimes-starttimes); % in seconds
        S_rippletimes.peakHeights = pkHeights;
        S_rippletimes.peakTimes = pkTimes; % in seconds
        S_rippletimes.peakWidths = pkWidths; 
        S_rippletimes.peakProms = pkProms;
        S_rippletimes.peakThr = peakThr;
        S_rippletimes.peakProm = peakProm;
        S_rippletimes.gaussWindowDur = gwinDur; % in ms
        S_rippletimes.minDipDur = minDipDur; % in ms
    
        % Store for each epoch
        rippletimes{1,e} = S_rippletimes;
    end
    % C_allripmeth{1,r} = riptimesEpochs;
    % To get rippletimes into the same format as Justin's rippletime
    % file, it must be embeded in an outer 1x1 cell array.
    rippletimes_MES = {};
    rippletimes_MES{1,1} = rippletimes;
    save(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/%s_direct/%srippletimes_MES.mat",rStr,short_name),"rippletimes_MES")
% end


% for m = 1:length(methodNames)
%     method = methodNames{m};
% 
%     if strcmp(method, "Justin")
%   
% 
% 
% 
% 
% 
% 
% 
% 
%     end
% end