% This script will dive a bit into ripples of both behavioral and rest
% epochs. 
% Note: Sometimes in this script I mistakenly say "multiunit activity" when I
% actually only mean cell activity. Multiunit activity refers to all the
% spikes detected on a tetrode, clustered to a cell or not, and is
% sometimes used for ripple detection. I am not doing this, and instead
% only using spikes from clustered cells.
% Michael Satchell Aug 2023

%% Load data

clearvars;

data_dir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
load_rats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','pos01'
filetypes = {'spikes01','tetinfo','pos01','sws01','rem01','rippletimes_MES','behavperform'};

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
            if iscell(file)
                C_alldata{ft,r} = file{1,1};
            else
                C_alldata{ft,r} = file;
            end
            fprintf("       Loaded file: %s  \n",  File_dir.name)
        end
    end
end

C_alldata = clip_17_epochs(C_alldata); % removes extra epoch data.


% ------------------------------------------
% ------------------------------------------
% ------------------------------------------
% List of all states to sort the spike data into. Note that adding
% ripple state adds quite a bit of calculation time. 
% WARNING: This does not change the order of how these states are stored in
% C_allstates or later FR_allStates. Thus, stateFiles MUST MATCH ORDER
% THESE FILES ARE LISTED IN filetypes.
stateFiles = {'sws','rem','ripple'};
stateNames = {'sws','rem','ripple','run','still'}; % Must match order states
% are added to C_allstates
% ------------------------------------------
% ------------------------------------------
% ------------------------------------------


sws_idx = find(contains(filetypes,'sws'));
if isempty(sws_idx)
    error("sws data must be loaded to run this analysis.")
end

rem_idx = find(contains(filetypes,'rem'));
if isempty(rem_idx)
    error("rem data must be loaded to run this analysis.")
end

spikes_idx = find(contains(filetypes,'spikes01')); 
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end

% cellinfo_idx = find(contains(filetypes,'cellinfo'));
% if isempty(cellinfo_idx)
%     error("cellinfo data must be loaded to run this analysis.")
% end

states_idx = find(contains(filetypes, stateFiles));
if isempty(states_idx)
    error("at least one of the following rest state data must be loaded: \n +" + ...
        "   sleep01, waking01, sws01, rem01, or rippletime01.")
end

% linf_idx = find(contains(filetypes,'linfields'));
% if isempty(linf_idx)
%     error("linfields data must be loaded to run this analysis.")
% end

pos_idx = find(contains(filetypes,'pos01'));
if isempty(pos_idx)
    error("pos01 data must be loaded to run this analysis.")
end

rt_idx = find(contains(filetypes,'rippletime'));
if isempty(rt_idx)
    error("rippletime data must be loaded to run this analysis.")
end

ti_idx = find(contains(filetypes,'tetinfo'));
if isempty(ti_idx)
    error("tetinfo data must be loaded to run this analysis.")
end

behper_idx = find(contains(filetypes,'behavperform'));
if isempty(behper_idx)
    error("behavperform data must be loaded to run this analysis.")
end

C_swsstate = C_alldata(sws_idx,:);
C_remstate = C_alldata(rem_idx,:);
C_allspikes = C_alldata(spikes_idx,:);
% C_allinfo = C_alldata(cellinfo_idx,:);
C_allstates = C_alldata(states_idx,:);
% C_alllinf = C_alldata(linf_idx,:);
C_runstate = create_runstate(C_alldata(pos_idx,:));
C_stillstate = create_stillstate(C_alldata(pos_idx,:));
C_allstates = [C_allstates; C_runstate; C_stillstate];
C_allriptimes = C_alldata(rt_idx,:);
C_alltetinfo = C_alldata(ti_idx,:);
C_allbehper = C_alldata(behper_idx,:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;


brainAreas = {'CA1','PFC'};

%% Calculates firing rates, state occurance times, and state durations
[FR_allStates,isPyrMat,isInhMat,Occ_allStates,Dur_allStates] = calc_meanrates(brainAreas,C_allstates,C_allspikes);


%% Exploratory ripple analyses

% % Plot total ripple duration across epochs
% ripDurs = zeros([size(C_allript,2), 17]);
% for r = 1:size(C_allript,2)
%     for e = 1:size(C_allript{1,r},2)
%         ripDurs(r,e) = C_allript{1,r}{1,e}.total_duration;
% 
%     end
% end
% 
% figure;
% hold on;
% plot(behEpochs,mean(ripDurs(:,behEpochs),1), Color=[0 0.4470 0.7410])
% plot(behEpochs,ripDurs(:,behEpochs)', '.', Color=[0 0.4470 0.7410])
% plot(restEpochs,mean(ripDurs(:,restEpochs),1), Color=[0.8500 0.3250 0.0980])
% plot(restEpochs,ripDurs(:,restEpochs)', '.', Color=[0.8500 0.3250 0.0980])
% % legend()
% title("Total Ripple Time per Epoch")
% ylabel("Epoch")
% xlabel("Total Duration (s)")




% % Plot LFP with ripple times in one rat/epoch. Also include SWS, still
% % and run times.
% % NOTE still states overlap with SWS states, so anytime there is a SWS
% % state that is gauranteed to also be a still state.
% r = 3; % rat number
% eStr = "01"; % epoch number
% e = str2num(eStr);
% riptetsStr = {'12','10','07','21','22','25'}; 
% riptets = [12,10,7, 21, 22,25];
% ydisp = [5,4,3,0,-3,-4];
% 
% 
% figure;
% hold on;
% ripple_labels = {};
% 
% for rt = 1:length(riptets)
%     % CHECK that rat number, epoch, and tetrode match the loaded file
%     load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/BG1_direct/EEG/BG1eegref01-%s-%s.mat",eStr,riptetsStr{rt}))
% 
%     % See main_explore on why I am using eegref instead of eeg.
% 
%     refData = eegref{1,1}{1,e}{1,riptets(rt)};
%     timeRange = (refData.starttime : 1/refData.samprate : refData.endtime)';
%     refV = refData.data/refData.voltage_scaling; % Remove scaling factor
%     refV = refV/1000; % Change from microvolts to millivolts.
% 
%     plot(timeRange, refV+ydisp(rt))
% 
% end
% 
% % Get ripple time data
% riptData = C_allriptimes{1,r}{1,e};
% for i = 1:length(riptData.starttime)
%     x_vertices = [riptData.starttime(i),riptData.endtime(i),riptData.endtime(i),riptData.starttime(i)];
%     y_vertices = [min(refV+ydisp(rt)),min(refV+ydisp(rt)),max(refV),max(refV)];
%     patch(x_vertices, y_vertices, [0.4940 0.1840 0.5560],'FaceAlpha', 0.6, 'EdgeColor','none')
% end
% 
% 
% % I need to plot velocity with the times I have designated as "run" to check that 
% % what I did worked.
% posData = C_alldata(find(contains(filetypes,'pos01')),:); % Grabs position file data.
% swsData = C_allstates(find(contains(stateNames,'sws')),:);
% stillData = C_allstates(find(contains(stateNames,'still')),:);
% runData = C_allstates(find(contains(stateNames,'run')),:);
% posStruct = posData{1,r}{1,e};
% swsStruct = swsData{1,r}{1,e};
% stillStruct = stillData{1,r}{1,e};
% runStruct = runData{1,r}{1,e};
% vel = posStruct.data(:,5);
% vel_sm = sgolayfilt(vel, 2, 31); % Smoothing the data.
% 
% % Plot velocity
% plot(posStruct.data(:,1), vel_sm)
% 
% % % Plot sws times
% % swsLabels = {};
% % for i = 1:length(swsStruct.starttime)
% %     x_vertices = [swsStruct.starttime(i),swsStruct.endtime(i),swsStruct.endtime(i),swsStruct.starttime(i)];
% %     y_vertices = [min(vel_sm),min(vel_sm),max(vel_sm)+1,max(vel_sm)+1];
% %     patch(x_vertices, y_vertices, [0 0.4470 0.7410],'FaceAlpha', 0.3, 'EdgeColor','none') % blue
% % 
% %     if i == 1
% %         swsLabels{i} = "sws times";
% %     else
% %         swsLabels{i} = "";
% %     end
% % end
% 
% % Plot still times
% stillLabels = {};
% for i = 1:length(stillStruct.starttime)
%     x_vertices = [stillStruct.starttime(i),stillStruct.endtime(i),stillStruct.endtime(i),stillStruct.starttime(i)];
%     y_vertices = [min(vel_sm),min(vel_sm),max(vel_sm),max(vel_sm)];
%     patch(x_vertices, y_vertices, [0 0.9 1],'FaceAlpha', 0.3, 'EdgeColor','none') % light blue
% 
%     if i == 1
%         stillLabels{i} = "still times";
%     else
%         stillLabels{i} = "";
%     end
% end
% 
% % Plot run times
% runLabels = {};
% for i = 1:length(runStruct.starttime)
%     x_vertices = [runStruct.starttime(i),runStruct.endtime(i),runStruct.endtime(i),runStruct.starttime(i)];
%     y_vertices = [min(vel_sm),min(vel_sm),max(vel_sm),max(vel_sm)];
%     patch(x_vertices, y_vertices, [0.6 0.9 0],'FaceAlpha', 0.3, 'EdgeColor','none') % light green
% 
%     if i == 1
%         runLabels{i} = "run times";
%     else
%         runLabels{i} = "";
%     end
% end
% 
% title(sprintf("Rat %d, Epoch %d",r,e))
% % ylabel("Voltage (mV)")
% xlabel("Time (s)")
% % legend(riptetsStr,swsLabels{:})




% % Plot LFP, ripple times and spikes of ripple tetrodes and all spikes on
% % those tetrodes. Also plot riple detection algorithm and peak finding.
% 
% rStr = "ZT2"; % rat name
% r = find(contains(load_rats,rStr));
% eStr = "03"; % epoch number
% e = str2num(eStr);
% 
% % All riptet ordering is based on epoch 3 for each animal:
% 
% % For ZT2 (might crash matlab)
% riptetsStr = {'10','11','12','13','14','15','16','17', '18','22','23', '19','20','21','24','25',...
%     '26','27','28','29','34','36', '31','32','33','35'};
% riptets = [10,11,12,13,14,15,16,17, 18,22,23, 19,20,21,24,25,...
%     26,27,28,29,34,36, 31,32,33,35];
% ydisp = [15:-1:8,  5:-1:3, 0:-1:-10 -13:-1:-16];
% 
% % % For ER1_NEW
% % % riptetsStr must be double-digit (i.e. must be '02', not '2').
% % riptetsStr = {'26','25','24','23', '11','07','21','08','15','13', '14'}; % tetrode to take LFP from, ideally displays most of the ripples
% % riptets = [26,25,24,23, 11,7,21,8,15,13, 14]; % ENSURE this matches riptetsStr
% % ydisp = [6,5,4,3, 0,-1,-2,-3,-4,-5, -8]; % y-displacements. Groupings of LFPs for clear plotting.
% 
% % % For KL8
% % riptetsStr = {'09','10', '25','24','15','21','23','22'}; 
% % riptets = [9,10, 25,24,15,21,23,22];
% % ydisp = [8,7, 4,3,2,1,0,-1];
% 
% % % For BG1
% % riptetsStr = {'12','10','07','21','22','25'}; 
% % riptets = [12,10,7, 21, 22,25];
% % ydisp = [5,4,3,0,-3,-4];
% 
% % % For JS14
% % riptetsStr = {'12','06','08','09','10','15','26'}; 
% % riptets = [12, 6,8,9,10,15,26];
% % ydisp = [7, 4,3,2,1,0,-1];
% 
% % % % For JS15
% % riptetsStr = {'05','06','07','08','09','11','21','24','25'};
% % riptets = [5,6,7,8,9,11,21,24,25];
% % ydisp = [5,4,3,2,1,0,-1,-2,-3];
% 
% % % For JS17
% % riptetsStr = {'06','05','24','22','11','12'};
% % riptets = [6,5,24, 22, 11,12];
% % ydisp = [6,5,4,1,-2,-3];
% 
% % % For JS21
% % riptetsStr = {'05','06','25','07', '08'};
% % riptets = [5,6,25,7, 8];
% % ydisp = [6,5,4,3,0];
% 
% % % For JS34
% % riptetsStr = {'06','07', '11','12', '22','24','25','09'};
% % riptets = [6,7, 11,12, 22,24,25,9];
% % ydisp = [8,7,4,3,0,-1,-2,-3];
% 
% 
% tetcolors = zeros(length(riptets),3);
% 
% 
% figure;
% hold on;
% ripple_labels = {};
% 
% for rt = 1:length(riptets)
%     % CHECK that rat number, epoch, and tetrode match the loaded file
%     load(sprintf("/media/msatchell/10TBSpinDisk/js_SingleDayExpt/%s_direct/EEG/%seegref01-%s-%s.mat",rStr,rStr,eStr,riptetsStr{rt}))
%     eegripple = load(sprintf("/media/msatchell/10TBSpinDisk/js_SingleDayExpt/%s_direct/EEG/%sripple01-%s-%s.mat",rStr,rStr,eStr,riptetsStr{rt}));
% 
%     % See main_explore on why I am using eegref instead of eeg.
% 
%     refData = eegref{1,1}{1,e}{1,riptets(rt)};
%     timeRange = (refData.starttime : 1/refData.samprate : refData.endtime)';
%     refV = refData.data/refData.voltage_scaling; % Remove scaling factor
%     refV = refV/1000; % Change from microvolts to millivolts.
% 
%     ripData = eegripple.ripple{1,1}{1,e}{1,riptets(rt)};
%     ripAmp = ripData.data(:,1)/ripData.voltage_scaling; % Ripple filtered LFP amplitude.
%     ripAmp = ripAmp/1000; % Change from microvolts to millivolts.
%     ripEnv = ripData.data(:,3); % Ripple filter envelope.
%     % Zscore the envelope of the ripple filter to determine significant
%     % ripples.
%     zEnv = zscore(double(ripEnv));
% 
%     if rt == 1 % define matrix to hold envelope traces. Envelope data should be the same length for all tetrodes.
%         riptetszEnvs = zeros([size(zEnv,1), length(riptets)]);
%     end
% 
%     riptetszEnvs(:,rt) = zEnv; % add envelope data
% 
%     tetline = plot(timeRange, refV+ydisp(rt));
%     c = get(tetline, 'Color');
%     tetcolors(rt,:) = c;
% 
% end
% 
% 
% % It is important to choose a good gaussian curve width for ripple
% % detection, because two things need to be considered: 1) cells in a ripple
% % do not fire simultaneously, so the width cannot be too small or no
% % overlap between the cell activity will be detected, and 2) cell bursting
% % may be a factor to consider, and by using a large width the difference
% % between a single spike and a small burst becomes small.
% % What I think is important for using multiunit activity in ripple analysis
% % is to incorporate the fact that ripples are special because they activate
% % a wide distribution of neurons. High cell activity and bursting can
% % happen outside of ripples, but the diversity in neuronal firing during the 
% % small window of a ripple is a somewhat signature trait, and should be
% % taken advantage of for ripple detection. To this effect, it is important
% % that for each neuron, gaussian curves do not sum when overlapping, as
% % this would increase the importance of bursting activity in the signal.
% % Instead, gaussians should be summed across different neurons so as to make the
% % signal large when there are many different neurons coincidentally active.
% % Unfortunately it is difficult for me to figure out how to convolve
% % without summing overlapping parts, so I will do the next closest thing
% % cut off the peaks at an amplitude of 1. So the overlapping tails of the gaussians
% % will sum, but the peaks will not.
% gwinlen = 151;
% gwin = gausswin(gwinlen); 
% % with a smaple rate of 1500 Hz, each time point is 2/3
% % of a milisecond apart, so to get a 4 ms gaussian for each spike, the
% % window must be 6 time points wide. Use an odd number so the gaussian is symmetric
% % about the spike. 30 points = 20ms. 150 = 100ms.
% 
% % NOTE I need to be careful when working with spike times. The end time for 
% % the LFP, from the one example I looked at, is 0.033 microseconds longer
% % than the end time of the range for the spike data. This means that the
% % timeRange arrays created using the samprate would misalign every 6 or 7
% % time points. So I can't simply interchange them. I will stick with the
% % LFP time range for now.
% 
% 
% % Plot neurons from all tetrodes
% nrnTets = []; % Tetrodes to plot neurons from.
% % Sort out CA1 tetrodes from PFC ones.
% for tet = 1:size(C_alltetinfo{1,r}{1,e},2)
%     tetinfo = C_alltetinfo{1,r}{1,e}{1,tet};
%     if isfield(tetinfo,"area") && strcmp(tetinfo.area,"CA1")
%         nrnTets(end+1) = tet;
%     end
% end
% 
% nrnCount = 0; % To count the number of plotted neurons in the raster for spacing purposes.
% nrnsGtrace = [];
% nrncolor = [0.8 0.8 0.8];
% for i = 1:length(nrnTets)
%     tet = nrnTets(i);
%     if ~isempty(C_allspikes{1,r}{1,e}{1,tet})
%         nrns = {C_allspikes{1,r}{1,e}{1,tet}{:}}; % cell with neuron structs.
%         if ~isempty(nrns)
%             idx = find(riptets == tet);
% 
%             for j = 1:size(nrns,2)
% 
%                 nrn = nrns{j};
% 
%                 if isfield(nrn,'meanrate')
%                     if nrn.meanrate <= 7
% 
%                         spikeTimes = nrn.data(:,1);
%                         % nrnTimeRange = nrn.timerange(1):1/refData.samprate:nrn.timerange(2);
% 
%                         yvals = zeros(size(spikeTimes)) + (min(ydisp)-2-0.5*nrnCount);
%                         plot(spikeTimes, yvals, Color=nrncolor, Marker="|", LineStyle="none")
% 
%                         % Convolve this window with the spike train for a neuron of
%                         % 1s and 0s. There was a problem that was occurring:
%                         % ismember needs exact equality between the timepoints in
%                         % spikeTmes and timeRange to identify matches, but the
%                         % spike times rarely line up with the times in timeRange,
%                         % regardless of whether the spike or LFP data start and end
%                         % times are used. So, ismembertol is used to fix this by
%                         % allowing a 1/samprate window of tolerance for determining
%                         % equality around each time point.
%                         isSpike = ismembertol(timeRange,spikeTimes, 1/(2*refData.samprate), DataScale=1);
%                         gtrace = conv(isSpike,gwin,"same");
%                         cutgtrace = gtrace > 1; % where to cut off summed peaks
%                         gtrace(cutgtrace) = 1; % reassign these values to 1.
%                         plot(timeRange,gtrace+(min(ydisp)-2-0.5*nrnCount),Color=nrncolor)
% 
%                         nrnCount = nrnCount + 1;
% 
%                         % add nrn gaussian trace to matrix
%                         nrnsGtrace(:,nrnCount) = gtrace;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% 
% 
% % Get ripple time data
% riptData = C_allriptimes{1,r}{1,e};
% for i = 1:length(riptData.starttime)
%     x_vertices = [riptData.starttime(i),riptData.endtime(i),riptData.endtime(i),riptData.starttime(i)];
%     y_vertices = [min(ydisp)-4-0.5*nrnCount,min(ydisp)-4-0.5*nrnCount,max(ydisp)+2,max(ydisp)+2];
%     patch(x_vertices, y_vertices, [0.4940 0.1840 0.5560],'FaceAlpha', 0.2, 'EdgeColor','none')
% end
% 
% % % My new rippletime data
% % riptData_new = riptimesEpochs{1,e};
% % for i = 1:length(riptData_new.starttime)
% %     x_vertices = [riptData_new.starttime(i),riptData_new.endtime(i),riptData_new.endtime(i),riptData_new.starttime(i)];
% %     y_vertices = [min(ydisp)-4-0.5*nrnCount,min(ydisp)-4-0.5*nrnCount,max(ydisp)+2,max(ydisp)+2];
% % patch(x_vertices, y_vertices, [0.2 1 0],'FaceAlpha', 0.2, 'EdgeColor','none')
% % end
% 
% % Plot run and still times as bars
% runStarttimes = C_runstate{1,r}{1,e}.starttime;
% runEndtimes = C_runstate{1,r}{1,e}.endtime;
% for i = 1:length(runStarttimes)
%     plot([runStarttimes(i), runEndtimes(i)], [min(ydisp)-8-0.5*nrnCount, min(ydisp)-8-0.5*nrnCount],...
%         LineWidth=5, Color=[0.6 0.9 0])
% end
% stillStarttimes = C_stillstate{1,r}{1,e}.starttime;
% stillEndtimes = C_stillstate{1,r}{1,e}.endtime;
% for i = 1:length(stillStarttimes)
%         plot([stillStarttimes(i), stillEndtimes(i)], [min(ydisp)-8.5-0.5*nrnCount, min(ydisp)-8.5-0.5*nrnCount],...
%         LineWidth=5, Color=[0 0.9 1])
% end
% 
% % Sum the gaussian traces from all neurons and zscore it.
% sumGtrace = sum(nrnsGtrace,2);
% zGtrace = zscore(sumGtrace);
% % plot(timeRange, (min(ydisp)-4-0.5*nrnCount)+zGtrace,  Color=[0.6 0.3 0]) % Brown
% 
% % Calculate zscore of sum of envelope zscores of all tetrodes and plot.
% % Each tetrode's z-scored envelope measures, relatively, the power of the ripple band 
% % for that tetrode, and if we add all the zEnvs together we get the
% % population relative zscore. Zscoring this sum shows the significance of
% % power changes for the all tetrodes while maintaining the difference in
% % mean ripple power for each tetrode.
% sumzEnv = sum(riptetszEnvs,2);
% popEnvzscore = zscore(sumzEnv); % Zscore of all tetrodes
% 
% zCum = zscore(zGtrace+popEnvzscore); % z-score the combined mulitunit and tetrode zscore traces.
% 
% peakThr = 6; % Number of standard deviations above the mean to detect peak
% peakSep = 0;  % Minimum time (in ms) between peaks.
% peakWidth = 10; % Minimum width (in ms) that a peak must have.
% peakProm = 6; % Minimum height (in stds) above neighboring signal a peak must have.
% 
% % Use the findpeaks function to find where the zscore passes above
% % threshold. Note that when the sample rate is provided, peak locations
% % and widths are returned in units of time.
% [pkHeights,pkTimes,pkWidths,pkProms] = findpeaks(zGtrace+popEnvzscore,timeRange, MinPeakHeight=peakThr,...
%     MinPeakDistance=peakSep/1000, MinPeakWidth=peakWidth/1000, MinPeakProminence=peakProm);
% % [heights,times,widths,proms] = findpeaks(zCom,timeRange, MinPeakHeight=peakThr,...
% %     MinPeakDistance=peakSep/1000, MinPeakWidth=peakWidth/1000, MinPeakProminence=peakProm);
% 
% plot(pkTimes, (min(ydisp)-4-0.5*nrnCount)+pkHeights, '.r',HandleVisibility="off")
% plot([min(timeRange),max(timeRange)],[(min(ydisp)-4-0.5*nrnCount)+peakThr,(min(ydisp)-4-0.5*nrnCount)+peakThr], '--b',HandleVisibility="off")
% plot([min(timeRange),max(timeRange)],[(min(ydisp)-4-0.5*nrnCount),(min(ydisp)-4-0.5*nrnCount)], '--k',HandleVisibility="off")
% plot([min(timeRange),max(timeRange)],[(min(ydisp)-4-0.5*nrnCount)+1,(min(ydisp)-4-0.5*nrnCount)+1], '--r',HandleVisibility="off")
% % plot(timeRange, min(ydisp)+popEnvzscore, 'k',HandleVisibility="off")
% plot(timeRange, (min(ydisp)-4-0.5*nrnCount)+zGtrace+popEnvzscore, 'k',HandleVisibility="off")
% % plot(timeRange, (min(ydisp)-4-0.5*nrnCount)+zCum, 'k',HandleVisibility="off")
% % plot(timeRange, (min(ydisp)-4-0.5*nrnCount)+popEnvzscore, Color=[1 0.4 0.3],HandleVisibility="off") % pink
% 
% title(sprintf("%s, Epoch %d | numRips=%d",rStr,e,length(riptData.starttime)))
% ylabel("Voltage (mV)")
% xlabel("Time (s)")




% % Plot voltage, envelope and zscores.
% eStr = "03"; % epoch number
% e = str2num(eStr);
% rStr = "JS21"; % rat name
% r = find(contains(load_rats,rStr));
% rt_idx = 1; % Index of riptet in C_riptets
% 
% load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/%s_direct/EEG/%seegref01-%s-%s.mat",rStr,rStr,eStr,C_riptets{1,r}{1,rt_idx}))
% load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/%s_direct/EEG/%sripple01-%s-%s.mat",rStr,rStr,eStr,C_riptets{1,r}{1,rt_idx}))
% 
% refData = eegref{1,1}{1,e}{1,str2num(C_riptets{1,r}{1,rt_idx})};
% timeRange = (refData.starttime : 1/refData.samprate : refData.endtime)';
% refV = refData.data/refData.voltage_scaling; % Remove scaling factor
% refV = refV/1000; % Change from microvolts to millivolts.
% 
% ripData = ripple{1,1}{1,e}{1,str2num(C_riptets{1,r}{1,rt_idx})};
% ripAmp = ripData.data(:,1)/ripData.voltage_scaling; % Ripple filtered LFP amplitude.
% ripAmp = ripAmp/1000; % Change from microvolts to millivolts.
% ripEnv = ripData.data(:,3); % Ripple filter envelope.
% 
% figure;
% hold on;
% % plot(timeRange, refV)
% % plot(timeRange, ripEnv)
% plot(timeRange, zscore(refV))
% plot(timeRange, zscore(double(ripEnv)))



%% Generate ripple info across animals

allRipRates = zeros(size(C_allriptimes,2),size(C_allriptimes{1,1},2)); % Ripple rates
allRipTotTime = zeros(size(C_allriptimes,2),size(C_allriptimes{1,1},2)); % Total ripple time normalized by the 
% length of SWS or still state time.
allRipAvgDur = zeros(size(C_allriptimes,2),size(C_allriptimes{1,1},2)); % average ripple duration
allEpochPer = zeros(size(C_allriptimes,2),size(C_allriptimes{1,1},2)); % Total performance across the entire epoch

for r = 1:size(C_allriptimes,2)
    % We want to calculate ripple rate as the number of ripples divided by
    % the total time available for ripples to occur, i.e. SWS and still
    % states. BUT still already encompasses all SWS times, so I should just
    % consider still times as periods during which ripples can happen,
    % otherwise I will be counting SWS times twice.
    % swsStates = C_swsstate{1,r};
    stillStates = C_stillstate{1,r};
    S_behper = C_allbehper{1,r}; % Performance data
    for e = 1:17
        
        S_stillstate = stillStates{1,e};
        stillDur = sum(S_stillstate.endtime - S_stillstate.starttime);
        
        if mod(e,2) == 0 % beh epochs
            % The trials are not split up by epoch, so I need to do it manually
            % from dayouttrials (which should be named epochouttrials).
            firstTrial = S_behper.dayouttrials(e/2,1);
            lastTrial = S_behper.dayouttrials(e/2,2);
            epochPer = S_behper.outreward(firstTrial:lastTrial); % Performance during this epoch
            allEpochPer(r,e) = sum(epochPer)/size(epochPer,1);
        end

        % fprintf("Still Dur: %.1f \n",sum(S_stillstate.endtime - S_stillstate.starttime))


        S_riptimes = C_allriptimes{1,r}{1,e};

        numRips = size(S_riptimes.starttime,1); % Number of independent ripple events
        allRipRates(r,e) = numRips/stillDur;

        allRipTotTime(r,e) = S_riptimes.total_duration/stillDur;

        allRipAvgDur(r,e) = S_riptimes.total_duration/numRips;
        
    end  

end

%% Plot ripple rate, duration, and cell firing rates




% % Plot average ripple rate, total ripple time, and ripple durations (with data points)
% figure
% hold on
% for r = 1:size(C_alldata,2)
%     plot(behEpochs, allRipRates(r,behEpochs), ".", Color='b', HandleVisibility="off")
%     plot(restEpochs, allRipRates(r,restEpochs), ".", Color='r', HandleVisibility="off")
% end
% plot(behEpochs, mean(allRipRates(:,behEpochs),1), Color='b',DisplayName="beh")
% plot(restEpochs, mean(allRipRates(:,restEpochs),1), Color='r',DisplayName="rest")
% title("Ripple Rates in Still and SWS")
% ylabel("Ripple Rate (ripples/sec)")
% xlabel("Epoch")
% legend()
% 
% figure
% hold on
% for r = 1:size(C_alldata,2)
%     plot(behEpochs, allRipTotTime(r,behEpochs), ".", Color="b", HandleVisibility="off")
%     plot(restEpochs, allRipTotTime(r,restEpochs), ".", Color='r', HandleVisibility="off")
% end
% plot(behEpochs, mean(allRipTotTime(:,behEpochs),1), Color='b',DisplayName="beh")
% plot(restEpochs, mean(allRipTotTime(:,restEpochs),1), Color='r',DisplayName="rest")
% title("% of Time Spent Rippling in Still and SWS")
% ylabel("(Total Ripple Time)/(SWS & Still Duration)")
% xlabel("Epoch")
% legend()
% 
% figure
% hold on
% for r = 1:size(C_alldata,2)
%     plot(behEpochs, allRipAvgDur(r,behEpochs), ".", Color="b", HandleVisibility="off")
%     plot(restEpochs, allRipAvgDur(r,restEpochs), ".", Color='r', HandleVisibility="off")
% end
% plot(behEpochs, mean(allRipAvgDur(:,behEpochs),1), Color='b',DisplayName="beh")
% plot(restEpochs, mean(allRipAvgDur(:,restEpochs),1), Color='r',DisplayName="rest")
% title("Avg Ripple Duration")
% ylabel("(Total Ripple Time (s))/(Num Ripples)")
% xlabel("Epoch")
% legend()


% % Plot individual rats ripple rate, duration, and total time rippling
% figure
% hold on
% for r = 1:size(C_alldata,2)
%     h1 = plot(behEpochs, allRipRates(r,behEpochs), Color='b', HandleVisibility="off");
%     h2 = plot(restEpochs, allRipRates(r,restEpochs), Color='r', HandleVisibility="off");
% end
% title("Ripple Rates in Still and SWS by Rat")
% ylabel("Ripple Rate (ripples/sec)")
% xlabel("Epoch")
% legend([h1,h2],{"beh","rest"})
% 
% figure
% hold on
% for r = 1:size(C_alldata,2)
%     h1 = plot(behEpochs, allRipTotTime(r,behEpochs), Color="b", HandleVisibility="off");
%     h2 = plot(restEpochs, allRipTotTime(r,restEpochs), Color='r', HandleVisibility="off");
% end
% title("% of Time Spent Rippling in Still and SWS by Rat")
% ylabel("(Total Ripple Time)/(SWS & Still Duration)")
% xlabel("Epoch")
% legend([h1,h2],{"beh","rest"})
% 
% figure
% hold on
% for r = 1:size(C_alldata,2)
%     h1 = plot(behEpochs, allRipAvgDur(r,behEpochs), Color="b", HandleVisibility="off");
%     h2 = plot(restEpochs, allRipAvgDur(r,restEpochs), Color='r', HandleVisibility="off");
% end
% title("Avg Ripple Duration by Rat")
% ylabel("(Total Ripple Time)/(Num Ripples)")
% xlabel("Epoch")
% legend([h1,h2],{"beh","rest"})


% % Plot the ripple information and performance on a seperate fig for each
% % rat
% for r = 1:size(C_alldata,2)
%     f1 = figure;
%     subplot(3,1,1)
%     sgtitle(sprintf("%s",load_rats{r}), interpreter="none")
%     hold on
%     h1 = plot(behEpochs, allRipRates(r,behEpochs), Color='b', HandleVisibility="off");
%     h2 = plot(restEpochs, allRipRates(r,restEpochs), Color='r', HandleVisibility="off");
%     h3 = plot(behEpochs, allEpochPer(r,behEpochs), Color='k', HandleVisibility="off");
%     title("Ripple Rates")
%     ylabel("Ripple Rate (ripples/sec)")
%     xlabel("Epoch")
%     legend([h1,h2,h3],{"beh","rest","perf"},Location="best")
% 
%     subplot(3,1,2)
%     hold on
%     h1 = plot(behEpochs, allRipTotTime(r,behEpochs), Color="b", HandleVisibility="off");
%     h2 = plot(restEpochs, allRipTotTime(r,restEpochs), Color='r', HandleVisibility="off");
%     % h3 = plot(behEpochs, allEpochPer(r,behEpochs), Color='k', HandleVisibility="off");
%     title("% of Time Spent Rippling in Still")
%     ylabel("(Total Ripple Time)/(SWS & Still Duration)")
%     xlabel("Epoch")
%     legend([h1,h2,h3],{"beh","rest"},Location="best")
% 
% 
%     subplot(3,1,3)
%     hold on
%     h1 = plot(behEpochs, allRipAvgDur(r,behEpochs), Color="b", HandleVisibility="off");
%     h2 = plot(restEpochs, allRipAvgDur(r,restEpochs), Color='r', HandleVisibility="off");
%     % h3 = plot(behEpochs, allEpochPer(r,behEpochs), Color='k', HandleVisibility="off");
%     title("Avg Ripple Duration")
%     ylabel("(Total Ripple Time)/(Num Ripples)")
%     xlabel("Epoch")
%     legend([h1,h2,h3],{"beh","rest"},Location="best")
% 
%     pause
%     close(f1)
% 
% end


% Plot firing rates for each state along side ripple measures


% % Plot rates by epoch and state, seperated into pyrimidal and inhibitory.
% % Pyramidal
% stateEpochF = figure;
% sgtitle("Mean Pyramidal FR (CA1 Ripples)")
% stateColors = [[0 0.4470 0.7410]; [1 0 0]; [0.5 0.1 1]; [0.6 0.9 0]; [0 0.9 1]]; % Color of data, will be plotted with stateNames order.
% for a = 1:size(FR_allStates,2)
% 
%     subplot(1,size(FR_allStates,2),a)
%     hold on;
% 
%     for s = 1:size(FR_allStates,1)
%         % Calculate mean firing rate, standard deviation, and SEM values for
%         % each rest epoch.
%         FRpyr = FR_allStates{s,a}; % get state FRs
%         FRpyr(~isPyrMat) = NaN; % Set non-pyramidal cells to NaN
%         FRmeans = mean(FRpyr,1,'omitnan');
%         FRstds = std(FRpyr,0,1,'omitnan');
%         numnrns = sum(~isnan(FRpyr),1);
%         FRsems = FRstds./sqrt(numnrns);
% 
%         FR_CIs = zeros(2,size(FRmeans,2));
%         CI_pct = 75; % Confidence interval value in percent.
%         for e = 1:size(FRmeans,2)
%             temp_CI = tinv([(1-CI_pct/100)/2, 1-(1-CI_pct/100)/2], numnrns(1,e)-1); 
%             FR_CIs(:,e) = bsxfun(@times, FRsems(e), temp_CI(:));
%         end
%         plotEpochs = 1:17;
%         plotEpochs = plotEpochs(~isnan(FRmeans));
%         errorbar(plotEpochs,FRmeans(plotEpochs),FR_CIs(1,plotEpochs),FR_CIs(2,plotEpochs),".", Color=stateColors(s,:))
%         plot(plotEpochs,FRmeans(plotEpochs), Color=stateColors(s,:), LineWidth=2, HandleVisibility='off')
%     end
% 
%     title(brainAreas{a})
%     ylabel("Mean Firing Rate (Hz)")
%     xlabel("Epoch")
%     L = legend([stateNames,"beh"]);
% 
% end
% % Inhibitory
% stateEpochF = figure;
% sgtitle("Mean Inhibitory FR (CA1 Ripples)")
% stateColors = [[0 0.4470 0.7410]; [1 0 0]; [0.5 0.1 1]; [0.6 0.9 0]; [0 0.9 1]]; % Color of data, will be plotted with stateNames order.
% for a = 1:size(FR_allStates,2)
% 
%     subplot(1,size(FR_allStates,2),a)
%     hold on;
% 
%     for s = 1:size(FR_allStates,1)
%         % Calculate mean firing rate, standard deviation, and SEM values for
%         % each rest epoch.
%         FRinh = FR_allStates{s,a}; % get state FRs
%         FRinh(~isInhMat) = NaN; % Set non-pyramidal cells to NaN
%         FRmeans = mean(FRinh,1,'omitnan');
%         FRstds = std(FRinh,0,1,'omitnan');
%         numnrns = sum(~isnan(FRinh),1);
%         FRsems = FRstds./sqrt(numnrns);
% 
%         FR_CIs = zeros(2,size(FRmeans,2));
%         CI_pct = 75; % Confidence interval value in percent.
%         for e = 1:size(FRmeans,2)
%             temp_CI = tinv([(1-CI_pct/100)/2, 1-(1-CI_pct/100)/2], numnrns(1,e)-1); 
%             FR_CIs(:,e) = bsxfun(@times, FRsems(e), temp_CI(:));
%         end
%         plotEpochs = 1:17;
%         plotEpochs = plotEpochs(~isnan(FRmeans));
%         errorbar(plotEpochs,FRmeans(plotEpochs),FR_CIs(1,plotEpochs),FR_CIs(2,plotEpochs),".", Color=stateColors(s,:))
%         plot(plotEpochs,FRmeans(plotEpochs), Color=stateColors(s,:), LineWidth=2, HandleVisibility='off')
%     end
% 
%     title(brainAreas{a})
%     ylabel("Mean Firing Rate (Hz)")
%     xlabel("Epoch")
%     L = legend([stateNames,"beh"]);
% 
% end


% Plot FR of beh or rest epochs seperately by state, seperated into pyrimidal and inhibitory.
for et = 1:2 % Loop through beh and rest 
    epochTypes = {'behavior','rest'};
    if strcmp(epochTypes{et}, 'behavior') % behavioral epochs
        plotEpochs = behEpochs;
        plotStateNames = {'ripple','run','still'};
        stateColors = [[0.5 0.1 1]; [0.6 0.9 0]; [0 0.9 1]];
        plotStateIdxs = find(contains(stateNames, plotStateNames));
    elseif strcmp(epochTypes{et}, 'rest') % Rest epochs
        plotEpochs = restEpochs;    
        plotStateNames = {'sws','rem','ripple','run','still'};
        stateColors = [[0 0.4470 0.7410]; [1 0 0]; [0.5 0.1 1]; [0.6 0.9 0]; [0 0.9 1]];
        plotStateIdxs = find(contains(stateNames, plotStateNames));
    end
    FR_States = FR_allStates(plotStateIdxs,:);

    PyrF = figure;
    sgtitle(sprintf("%s Mean Pyramidal FR (CA1 Ripples)",epochTypes{et}))
    for a = 1:size(FR_States,2)
       
        subplot(1,size(FR_States,2),a)
        hold on;
        for s = 1:size(FR_States,1)
            % Calculate mean firing rate, standard deviation, and SEM values for
            % each rest epoch.
            FRpyr = FR_States{s,a}; % get state FRs
            FRpyr(~isPyrMat) = NaN; % Set non-pyramidal cells to NaN
            FRmeans = mean(FRpyr,1,'omitnan');
            FRstds = std(FRpyr,0,1,'omitnan');
            numnrns = sum(~isnan(FRpyr),1);
            FRsems = FRstds./sqrt(numnrns);
            
            FR_CIs = zeros(2,size(FRmeans,2));
            CI_pct = 75; % Confidence interval value in percent.
            for e = 1:size(FRmeans,2)
                temp_CI = tinv([(1-CI_pct/100)/2, 1-(1-CI_pct/100)/2], numnrns(1,e)-1); 
                FR_CIs(:,e) = bsxfun(@times, FRsems(e), temp_CI(:));
            end
            errorbar(plotEpochs,FRmeans(plotEpochs),FR_CIs(1,plotEpochs),FR_CIs(2,plotEpochs),".", Color=stateColors(s,:))
            plot(plotEpochs,FRmeans(plotEpochs), Color=stateColors(s,:), LineWidth=2, HandleVisibility='off')
        end
    
        title(brainAreas{a})
        ylabel("Mean Firing Rate (Hz)")
        xlabel("Epoch")
        L = legend(plotStateNames);
    
    end
    % Inhibitory
    InhF = figure;
    sgtitle(sprintf("%s Mean Inhibitory FR (CA1 Ripples)",epochTypes{et}))
    for a = 1:size(FR_States,2)
       
        subplot(1,size(FR_States,2),a)
        hold on;
        for s = 1:size(FR_States,1)
            % Calculate mean firing rate, standard deviation, and SEM values for
            % each rest epoch.
            FRinh = FR_States{s,a}; % get state FRs
            FRinh(~isInhMat) = NaN; % Set non-pyramidal cells to NaN
            FRmeans = mean(FRinh,1,'omitnan');
            FRstds = std(FRinh,0,1,'omitnan');
            numnrns = sum(~isnan(FRinh),1);
            FRsems = FRstds./sqrt(numnrns);
            
            FR_CIs = zeros(2,size(FRmeans,2));
            CI_pct = 75; % Confidence interval value in percent.
            for e = 1:size(FRmeans,2)
                temp_CI = tinv([(1-CI_pct/100)/2, 1-(1-CI_pct/100)/2], numnrns(1,e)-1); 
                FR_CIs(:,e) = bsxfun(@times, FRsems(e), temp_CI(:));
            end
            % plotEpochs = 1:17;
            % plotEpochs = plotEpochs(~isnan(FRmeans));
            errorbar(plotEpochs,FRmeans(plotEpochs),FR_CIs(1,plotEpochs),FR_CIs(2,plotEpochs),".", Color=stateColors(s,:))
            plot(plotEpochs,FRmeans(plotEpochs), Color=stateColors(s,:), LineWidth=2, HandleVisibility='off')
        end
    
        title(brainAreas{a})
        ylabel("Mean Firing Rate (Hz)")
        xlabel("Epoch")
        L = legend(plotStateNames);
    
    end


end



% PLOT FRACTION OF CELLS PARTICIPATING IN RIPPLES! 





