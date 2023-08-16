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


% % Plot LFP with ripple times above a raster plot of all cells in one rat/epoch
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


% I should also order the place cells based on peak firing location along
% one (or all) of the trajectories.
r = 2; % rat number
eStr = "06"; % epoch number
riptetsStr = {'26','25','24','23', '11','07','21','08','15','13', '14'}; % tetrode to take LFP from, ideally displays most of the ripples
e = str2num(eStr);
riptets = [26,25,24,23, 11,7,21,8,15,13, 14]; % ENSURE this matches riptetsStr
ydisp = [6,5,4,3, 0,-1,-2,-3,-4,-5, -8]; % y-displacements. Groupings of LFPs for clear plotting.

tetcolors = zeros(length(riptets),3);

figure;
hold on;
ripple_labels = {};

for rt = 1:length(riptets)
    % CHECK that rat number, epoch, and tetrode match the loaded file
    load(sprintf("/mnt/10TBSpinDisk/js_SingleDayExpt/ER1_NEW_direct/EEG/ER1eegref01-%s-%s.mat",eStr,riptetsStr{rt}))
    
    % See main_explore on why I am using eegref instead of eeg.
    
    refData = eegref{1,1}{1,e}{1,riptets(rt)};
    timeRange = (refData.starttime : 1/refData.samprate : refData.endtime)';
    refV = refData.data/refData.voltage_scaling; % Remove scaling factor
    refV = refV/1000; % Change from microvolts to millivolts.
    
    tetline = plot(timeRange, refV+ydisp(rt));
    c = get(tetline, 'Color');
    tetcolors(rt,:) = c;

end

nrnTets = [26,25,24,23, 11,7,21,8,15,13, 14]; % Tetrodes to plot neurons from.
nrnCount = 0; % To count the number of plotted neurons in the raster for spacing purposes.
for i = 1:length(nrnTets)
    tet = nrnTets(i);
    nrns = [C_allspikes{1,r}{1,e}{1,tet}{:}]; % struct with rows of neurons.
    if ~isempty(nrns)
        meanRates = [nrns.meanrate]';
        isI = meanRates > 7; % Inhibitory neurons.
        nrns = nrns(~isI); % Remove inhibitory neurons from struct.

        idx = find(riptets == tet);
    
        for j = 1:length(nrns)
            nrn = nrns(j);
            spikeTimes = nrn.data(:,1);
            
            yvals = zeros(size(spikeTimes)) + (min(ydisp)-2-0.5*nrnCount);
            plot(spikeTimes, yvals, Color=tetcolors(idx,:), Marker="|", LineStyle="none")

            nrnCount = nrnCount + 1;
        end
    end
end

% Get ripple time data
riptData = C_allript{1,r}{1,e};
for i = 1:length(riptData.starttime)
    x_vertices = [riptData.starttime(i),riptData.endtime(i),riptData.endtime(i),riptData.starttime(i)];
    y_vertices = [min(ydisp)-2-0.5*nrnCount,min(ydisp)-2-0.5*nrnCount,max(refV),max(refV)];
    patch(x_vertices, y_vertices, [0.4940 0.1840 0.5560],'FaceAlpha', 0.2, 'EdgeColor','none')
end
title(sprintf("Rat %d, Epoch %d",r,e))
ylabel("Voltage (mV)")
xlabel("Time (s)")
legend(riptetsStr)

