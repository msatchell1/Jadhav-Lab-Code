clearvars;

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
load_rats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01'.
filetypes = {'cellinfo','spikes01','pos01','sleep01','waking01','sws01','rem01'};

C_alldata = {}; % Cell array to hold data for all rats. If multiple filetypes 
% are loaded, each row holds a different file type, ordered in the same
% order as the elements in filetypes.

disp("Loading new animal data... ")
for r = 1:length(load_rats)
    fprintf("Loaded animal: %s \n",load_rats{r})

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




%% Split the firing rates into the following states during rest epochs:
% Behavior(on track), wake during sleep sessions, sleep, NREM, and REM.

spikes_idx = find(contains(filetypes,'spikes01')); % Gets row index of spike data in C_alldata.
% Note that indexing like this is only supposed to give one index as a
% result.
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end

cellinfo_idx = find(contains(filetypes,'cellinfo'));
if isempty(cellinfo_idx)
    error("cellinfo data must be loaded to run this analysis.")
end

% List of all states to sort the spike data into. Note that adding
% rippletime01 adds quite a bit of calculation time. 
% WARNING: This does not change the order of how these states are stored in
% C_allstates or later FR_allStates. Thus, stateNames MUST MATCH ORDER
% THESE FILES ARE LISTED IN filetypes.
stateNames = {'pos','sleep','waking','sws','rem'};

states_idx = find(contains(filetypes, stateNames));

% for s = 1:length(stateNames)
%     states_idx(s) = find(contains(filetypes, stateNames{s}));
% end
if isempty(states_idx)
    error("at least one of the following rest state data must be loaded: \n +" + ...
        "   sleep01, waking01, sws01, rem01, or rippletime01.")
end

brainAreas = {'CA1','PFC'}; % Brain areas to split data into. Should be only CA1 and PFC

% Extract spike data
C_allspikes = C_alldata(spikes_idx,:);

% Cell info data so I can determine brain region of each nrn. I could do a
% separate loop here and add an 'area' field to the neurons in C_allspikes.
% I could also just treat C_allinfo like C_allspikes and flatten the cells
% across tetrodes into one array inside the main for-loop, then grab the
% area from C_allinfo. For now I am going to choose the first option.
C_allinfo = C_alldata(cellinfo_idx,:);

% Add 'area' field to the spike data for later sorting.
for r = 1:length(C_allspikes)
    for e = 1:length(C_allspikes{1,r})
        for tet = 1:length(C_allspikes{1,r}{1,e})
            if ~isempty(C_allspikes{1,r}{1,e}{1,tet})
                for nrn = 1:length(C_allinfo{1,r}{1,e}{1,tet})
                    if ~isempty(C_allspikes{1,r}{1,e}{1,tet}{1,nrn})

                        % If the cell exists in the spike file it should
                        % exist in the cellinfo file.
                        nrnCellinfo = C_allinfo{1,r}{1,e}{1,tet}{1,nrn};
                        if isfield(nrnCellinfo,'area')
                            C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.area = nrnCellinfo.area; % Create field in spike struct.
                           
                        end
                    end
                end
            end
        end
    end
end


% Store all state data in one array. The order of the data
% stored in rows is that of stateNames.
C_allstates = C_alldata(states_idx,:);
% Identify times at which velocity is above 4 cm/s and fill a cell array
% full of structs that matches that of the other states like REM, NREM,
% etc. This is done by create_runstate.m.
pos_idx = find(contains(filetypes,'pos01'));
C_runstate = create_runstate(C_alldata(pos_idx,:));
% Now overwrite the pos data n C_allstates with the run state information
% in C_runstates.
for r = 1:size(C_runstate,2)
    C_allstates{find(contains(stateNames,'pos')),r} = C_runstate{1,r};
end

% Cell array to hold matrices for each brain region (in brainArea order).
% Matrices contain mean spike rates by epoch, not discriminating by state. 
% Rows are neuron identity, columns are epochs. This is for all rats.
FR_byEpoch = cell(1,length(brainAreas));

% Calculate and store mean firing rate across the entire epoch for each
% epoch, rat, and brain area.
for r = 1:length(C_allspikes) 
    
    for a = 1:length(brainAreas) % I could eliminate this loop by using strfind(brainAreas) to index...

        % Creates a matrix to hold spike rate data for one rat. Dims (num nrns)
        % x (num epochs). The first epoch is used to find the number of nrns
        % because this does not change across epochs as stated above. Third
        % dimension is area in order brainAreas.
        ratRates = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));

        for e = 1:length(C_allspikes{1,r})
    
            nrnsAlltets = [C_allspikes{1,r}{1,e}{:}]; % Combining nrn data from all tets.
    
            for nrn = 1:length(nrnsAlltets)
    
                % isfield returns false if the struct does not exist
                if isfield(nrnsAlltets{nrn},'data')
                    
                    if strcmp(nrnsAlltets{nrn}.area, brainAreas{a})
                        % I checked and confirmed that calculating the mean firing rate 
                        % this way matches the .meanrate field value for each nrn (only off by 0.0001 for some nrns).
                        ratRates(nrn,e) = size(nrnsAlltets{nrn}.data,1)/(nrnsAlltets{nrn}.timerange(2)-nrnsAlltets{nrn}.timerange(1));
                    end
                end
            end
        end
        FR_byEpoch{a} = [FR_byEpoch{a}; ratRates]; 
    end
end

% Separating into rest and behavior epoch data.
behEpochs = 2:2:size(FR_byEpoch{1},2);
restEpochs = 1:2:size(FR_byEpoch{1},2);
FRbeh = {}; FRrest = {};
for a = 1:size(FR_byEpoch,2)
    FRbeh{a} = FR_byEpoch{a}(:,behEpochs);
    FRrest{a} = FR_byEpoch{a}(:,restEpochs);
end


% Gets epochs that have data for all states. For most states this should be all odd
% epochs for the data I am sorting here, because all the states
% happen during the rest epochs except for rippletimes. I will detect non-empty
% epochs for generality. Note this assumes all rats have the same number of
% epochs as rat 1, which should be true when clip_17_epochs.m has been
% used on C_allstates.
nonemptyEpochs = zeros(size(C_allstates,1),length(C_allstates{1,1}));
for s = 1:size(C_allstates,1)
    nonemptyEpochs(s,:) = ~cellfun(@isempty, C_allstates{s,1});
end





% To hold firing rates for all states held in C_allstates
FR_allStates = cell(size(C_allstates,1),length(brainAreas));
% Time duration of states by epoch. (num states) x  (num rats) x (num epochs).
Dur_allStates = zeros(size(C_allstates,1),size(C_allstates,2),size(nonemptyEpochs,2));
Dur_ratRef = zeros(size(C_allstates,1),size(C_allstates,2),size(nonemptyEpochs,2)); % Array to reference rat number
% All the start and end times of all occurances. (num states) x (num rats).
% Then within each cell are cells for each epoch, and within those is a
% (num occurances) x 2 matrix. Col 1 is starttimes, col 2 endtimes.
Occ_allStates = cell(size(C_allstates,1),size(C_allstates,2));

% Main loop to calculate mean firing rate for all states in all brain areas.
for a = 1:length(brainAreas)

    for s = 1:size(C_allstates,1)
    
        % Gets epochs that have data for this state. This should be identical across rats.
        nonempty = zeros(1,length(C_allstates{1,1}));
    
        for r = 1:size(C_allstates,2) 
            % Indicates which epochs contain data
            nonempty(1,:) = ~cellfun(@isempty, C_allstates{s,r});
            % To hold mean firing rate data across all epochs for a single rat.
            ratMeans = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));
            % To hold state occurances for all epochs.
            occEpochs = cell(1,length(C_allstates{s,r}));
    
            for e = 1:length(C_allstates{s,r})
                
                if nonempty(1,e) == 1
                % Now the real deal. For every non-empty epoch, each state type (rem, rippletime,
                % sleep, sws, and waking) will loop through the times
                % that that state occurs, counting the number of spikes
                % that occur during the duration of that occurance. In the end
                % these spike counts will be summed up and divided by the
                % total duration of that state to get the mean rate.
                
                epochData = C_allstates{s,r}{1,e}; % Data on a single rat, state, and epoch.
                nrnsAlltets = [C_allspikes{1,r}{1,e}{:}]; % Combining nrn spike data from all tets.
                STs_inState = cell(1,length(nrnsAlltets)); % Cell of spike times that occur within the given state.
                % Columns are neurons.
    
                for o = 1:size(epochData.starttime,1) % Loop through occurances of that state.
    
                    for nrn = 1:length(nrnsAlltets)
    
                        % isfield returns false if the struct does not exist
                        if isfield(nrnsAlltets{nrn},'data') && ~isempty(nrnsAlltets{nrn}.data)...
                                && strcmp(nrnsAlltets{nrn}.area, brainAreas{a})
            
                            STs = nrnsAlltets{nrn}.data(:,1); % Time of all spikes for that nrn.
    
                            % Logical 1s and 0s determining if spikes fall within the state occurance window.
                            inOcc = (epochData.starttime(o) <= STs) & (STs <= epochData.endtime(o));
                            % Appends data from occurance to the state.
                            STs_inState{1,nrn} = [STs_inState{1,nrn}; STs(inOcc)];
        
                        end
                    end
                end
                    
                % Now that spikes have been collected for each nrn for all
                % occurances of the state, the mean rates for that state are calculated.
                % This means dividing the number of spikes by the total
                % duration of the state.
                occTimes = [epochData.starttime, epochData.endtime];
                occDurs = epochData.endtime - epochData.starttime;
                stateDur = sum(occDurs);  % Total time spent in state.
                Dur_allStates(s,r,e) = stateDur;
                Dur_ratRef(s,r,e) = r;
                FRmeans = NaN(size(STs_inState)); % To hold mean firing rate of nrns.
                for nrn = 1:size(STs_inState,2)
                    if ~isempty(STs_inState{1,nrn}) % Important to leave nrns with no spikes as NaNs.
                        meanFR = length(STs_inState{1,nrn})/stateDur; % Calculate mean firing rate
                        FRmeans(1,nrn) =  meanFR;
                    end
                end
                 
                ratMeans(:,e) = FRmeans; % Stores all FR means for that epoch.
                occEpochs{1,e} = occTimes;
               
                end
            end

            FR_allStates{s,a} = [FR_allStates{s,a}; ratMeans];    
            Occ_allStates{s,r} = occEpochs;
        end
    end
end


stateNames{find(contains(stateNames,'pos'))} = 'run'; % Change name to run 
% because that is a better description of the state.


%% I need to plot velocity with the times I have designated as "run" to check that 
% what I did worked.



%% Begin FR scatter plot analysis

% logical to plot all individual plots.
plotindv = 1;

stateColors = [[0.6 0.9 0]; [0 0.4470 0.7410]; [1 0.8 0]; [0.5 0.1 1]; [1 0 0]];

% SWS vs. all states
for a = 1:length(brainAreas)
    figure;
    sgtitle(sprintf("Neuron FRs | SWS vs. All States | %s",brainAreas{a}))
    for s = 1:size(FR_allStates,1)
        
        % Form vectors of firing rates for the two states being plotted against
        % each other, making sure to perserve the cell identity in each vector.
        % Each vector should hold info from all epochs.
        stateidx = find(contains(stateNames,'sws'));
    
        
        swsData = FR_allStates{stateidx,a};
        vsData = FR_allStates{s,a}; % Data of state to be plotted against.
    
        swsData = swsData(:,restEpochs);
        vsData = vsData(:,restEpochs);
        % NOTE I need to leave the NaNs in these matrices to preserve the
        % neuron IDs. Some neurons have NaN in one matrix and not the
        % other, so removing NaNs messes up the matching of x and y values
        % between the two lists. 
        
        subplot(ceil(size(FR_allStates,1)/2),2,s);
        hold on;

        numg = 1; % Number of groups to evenly split the data into.
        groups = randi([1,numg],[length(swsData(:)),1]);
        clrs = [];
        if numg > 1
            for g = 1:numg-1
                clrs = [clrs; [0.7 0.7 0.7]]; % Assign all but last group to color grey
            end
        end
        clrs = [clrs; stateColors(s,:)]; 

        gscatter(gca, swsData(:),vsData(:),randinds, clrs, '.', 1, doleg='off');
        plot(min(vsData(:)):max(vsData(:)), min(vsData(:)):max(vsData(:)),'k')
        xlabel("SWS Mean FR (Hz)")
        ylabel(sprintf("%s Mean FR (Hz)", stateNames{s}))
        set(gca, 'XScale', 'log', 'YScale', 'log');
    end
end


if plotindv

    







