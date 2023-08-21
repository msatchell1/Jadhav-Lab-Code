% This script will dive a bit into ripples of both behavioral and rest
% epochs. 

%% Load data

clearvars;

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
load_rats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01'.
filetypes = {'rippletime','linfields01','cellinfo','spikes01','pos01','sws01','rem01'};

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


stateFiles = {'sws','rem'};

stateNames = {'sws','rem','run','still'};

brainAreas = {'CA1','PFC'}; 

spikes_idx = find(contains(filetypes,'spikes01')); 
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end

cellinfo_idx = find(contains(filetypes,'cellinfo'));
if isempty(cellinfo_idx)
    error("cellinfo data must be loaded to run this analysis.")
end

states_idx = find(contains(filetypes, stateFiles));
if isempty(states_idx)
    error("at least one of the following rest state data must be loaded: \n +" + ...
        "   sleep01, waking01, sws01, rem01, or rippletime01.")
end

linf_idx = find(contains(filetypes,'linfields'));
if isempty(linf_idx)
    error("linfields data must be loaded to run this analysis.")
end

pos_idx = find(contains(filetypes,'pos01'));
if isempty(pos_idx)
    error("pos01 data must be loaded to run this analysis.")
end

rt_idx = find(contains(filetypes,'rippletime'));
if isempty(rt_idx)
    error("rippletime data must be loaded to run this analysis.")
end

C_allspikes = C_alldata(spikes_idx,:);
C_allinfo = C_alldata(cellinfo_idx,:);
C_allstates = C_alldata(states_idx,:);
C_alllinf = C_alldata(linf_idx,:);
C_runstate = create_runstate(C_alldata(pos_idx,:));
C_stillstate = create_stillstate(C_alldata(pos_idx,:));
C_allstates = [C_allstates; C_runstate; C_stillstate];
C_allript = C_alldata(rt_idx,:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;




%% Basic ripple analyses

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


% % Plot LFP with ripple times in one rat/epoch
% r = 2; % rat number
% eStr = "03"; % epoch number
% riptetsStr = {'26','25','24','23', '11','07','21','08','15','13', '14'}; % tetrode to take LFP from, ideally displays most of the ripples
% e = str2num(eStr);
% riptets = [26,25,24,23, 11,7,21,8,15,13, 14]; % ENSURE this matches riptetsStr
% ydisp = [6,5,4,3, 0,-1,-2,-3,-4,-5, -8]; % Groupings of LFPs for clear plotting.
% 
% figure;
% hold on;
% ripple_labels = {};
% 
% for rt = 1:length(riptets)
%     % CHECK that rat number, epoch, and tetrode match the loaded file
%     load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/ER1_NEW_direct/EEG/ER1eegref01-%s-%s.mat",eStr,riptetsStr{rt}))
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
% riptData = C_allript{1,r}{1,e};
% for i = 1:length(riptData.starttime)
%     x_vertices = [riptData.starttime(i),riptData.endtime(i),riptData.endtime(i),riptData.starttime(i)];
%     y_vertices = [min(refV+ydisp(rt)),min(refV+ydisp(rt)),max(refV),max(refV)];
%     patch(x_vertices, y_vertices, [0.4940 0.1840 0.5560],'FaceAlpha', 0.6, 'EdgeColor','none')
% end
% title(sprintf("Rat %d, Epoch %d",r,e))
% ylabel("Voltage (mV)")
% xlabel("Time (s)")
% legend(riptetsStr)


% % Plot LFP, ripple times and spikes of ripple tetrodes and all spikes on
% % those tetrodes.
% rStr = "ZT2"; % rat name
% r = find(contains(load_rats,rStr));
% eStr = "03"; % epoch number
% e = str2num(eStr);
% 
% % All riptet ordering is based on epoch 3 for each animal:
% 
% % % For ZT2
% % riptetsStr = {'10','11','12','13','14','15','16','17', '18','22','23', '19','20','21','24','25',...
% %     '26','27','28','29','34','36', '31','32','33','35'};
% % riptets = [10,11,12,13,14,15,16,17, 18,22,23, 19,20,21,24,25,...
% %     26,27,28,29,34,36, 31,32,33,35];
% % ydisp = [15:-1:8,  5:-1:3, 0:-1:-10 -13:-1:-16];
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
% % % For JS15
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
% figure;
% hold on;
% ripple_labels = {};
% 
% for rt = 1:length(riptets)
%     % CHECK that rat number, epoch, and tetrode match the loaded file
%     load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/%s_direct/EEG/%seegref01-%s-%s.mat",rStr,rStr,eStr,riptetsStr{rt}))
% 
%     % See main_explore on why I am using eegref instead of eeg.
% 
%     refData = eegref{1,1}{1,e}{1,riptets(rt)};
%     timeRange = (refData.starttime : 1/refData.samprate : refData.endtime)';
%     refV = refData.data/refData.voltage_scaling; % Remove scaling factor
%     refV = refV/1000; % Change from microvolts to millivolts.
% 
%     tetline = plot(timeRange, refV+ydisp(rt));
%     c = get(tetline, 'Color');
%     tetcolors(rt,:) = c;
% 
% end
% 
% nrnTets = riptets; % Tetrodes to plot neurons from.
% nrnCount = 0; % To count the number of plotted neurons in the raster for spacing purposes.
% for i = 1:length(nrnTets)
%     tet = nrnTets(i);
%     if ~isempty(C_allspikes{1,r}{1,e}{1,tet})
%         nrns = [C_allspikes{1,r}{1,e}{1,tet}{:}]; % struct with rows of neurons.
%         if ~isempty(nrns)
%             meanRates = [nrns.meanrate]';
%             isI = meanRates > 7; % Inhibitory neurons.
%             nrns = nrns(~isI); % Remove inhibitory neurons from struct.
% 
%             idx = find(riptets == tet);
% 
%             for j = 1:length(nrns)
%                 nrn = nrns(j);
%                 spikeTimes = nrn.data(:,1);
% 
%                 yvals = zeros(size(spikeTimes)) + (min(ydisp)-2-0.5*nrnCount);
%                 plot(spikeTimes, yvals, Color=tetcolors(idx,:), Marker="|", LineStyle="none")
% 
%                 nrnCount = nrnCount + 1;
%             end
%         end
%     end
% end
% 
% % Get ripple time data
% riptData = C_allript{1,r}{1,e};
% for i = 1:length(riptData.starttime)
%     x_vertices = [riptData.starttime(i),riptData.endtime(i),riptData.endtime(i),riptData.starttime(i)];
%     y_vertices = [min(ydisp)-2-0.5*nrnCount,min(ydisp)-2-0.5*nrnCount,max(ydisp)+2,max(ydisp)+2];
%     patch(x_vertices, y_vertices, [0.4940 0.1840 0.5560],'FaceAlpha', 0.2, 'EdgeColor','none')
% end
% title(sprintf("%s, Epoch %d",rStr,e))
% ylabel("Voltage (mV)")
% xlabel("Time (s)")
% legend(riptetsStr)





%% Create new rippletime files with improved SWR detection

% First, I need to come up with a ripple detection algorithm. For doing
% this, I think for now I will use the same riptets that Justin used for
% determining his rippletimes.


% A cell array that holds riptets in their plotting order for all rats.
% NOTE: Assumed order of rats in 'load_rats' is: 
% {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'}
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


eStr = "03"; % epoch number
e = str2num(eStr);
rStr = "JS21"; % rat name
r = find(contains(load_rats,rStr));
rt_idx = 1; % Index of riptet in C_riptets

load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/%s_direct/EEG/%seegref01-%s-%s.mat",rStr,rStr,eStr,C_riptets{1,r}{1,rt_idx}))
load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/%s_direct/EEG/%sripple01-%s-%s.mat",rStr,rStr,eStr,C_riptets{1,r}{1,rt_idx}))

refData = eegref{1,1}{1,e}{1,str2num(C_riptets{1,r}{1,rt_idx})};
timeRange = (refData.starttime : 1/refData.samprate : refData.endtime)';
refV = refData.data/refData.voltage_scaling; % Remove scaling factor
refV = refV/1000; % Change from microvolts to millivolts.

ripData = ripple{1,1}{1,e}{1,str2num(C_riptets{1,r}{1,rt_idx})};
ripAmp = ripData.data(:,1)/ripData.voltage_scaling; % Ripple filtered LFP amplitude.
ripAmp = ripAmp/1000; % Change from microvolts to millivolts.
ripEnv = ripData.data(:,3); % Ripple filter envelope.



figure;
hold on;
% plot(timeRange, refV)
% plot(timeRange, ripEnv)
plot(timeRange, zscore(refV))
plot(timeRange, zscore(double(ripEnv)))

% for r = 1:length(load_rats)
% 
%     rStr = load_rats{1,r};
% 
% 
% load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/%s_direct/EEG/%seegref01-%s-%s.mat",rStr,rStr,eStr,riptetsStr{rt}))
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% end


% put a lowpass filter on the LFP data to try and remove the spikes 