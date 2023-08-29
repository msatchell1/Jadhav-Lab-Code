clearvars;

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
load_rats = {'KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','rippletimes_MES.
filetypes = {'linfields01','cellinfo','spikes01','pos01','sws01','rem01','rippletimes_MES'};

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




%% Split the firing rates into the following states during rest epochs:
% Behavior(on track), wake during sleep sessions, sleep, NREM, and REM.

% ----------------------------------
% To remove states from the analysis, simply take them out of stateFiles
% and stateNames. They do not need to be removed from filetypes.

% List of all states to sort the spike data into. Note that adding
% rippletime01 adds quite a bit of calculation time. 
% WARNING: This does not change the order of how these states are stored in
% C_allstates or later FR_allStates. Thus, stateNames MUST MATCH ORDER
% THESE FILES ARE LISTED IN filetypes.
% stateFiles = {'sleep','waking','sws','rem'};
stateFiles = {'sws','rem','ripple'};

% Names to be used when plotting of the different states. Note, this array
% can be larger than stateFiles because some stateFiles can be used to
% produce multiple states, such as pos01 being used for 'run' and 'still'.
% This list must match the order of states in C_allstates.
% stateNames = {'sleep','wake','sws','rem','run','still'};
stateNames = {'sws','rem','ripple','run','still'};

brainAreas = {'CA1','PFC'}; % Brain areas to split data into. Should be only CA1 and PFC

% -----------------------------------------------------

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

states_idx = find(contains(filetypes, stateFiles));
if isempty(states_idx)
    error("at least one of the following rest state data must be loaded: \n +" + ...
        "   sleep01, waking01, sws01, rem01, or rippletime01.")
end

linf_idx = find(contains(filetypes,'linfields'));
if isempty(linf_idx)
    error("linfields data must be loaded to run this analysis.")
end


% Extract spike data
C_allspikes = C_alldata(spikes_idx,:);

% Cell info data so I can determine brain region of each nrn. I could do a
% separate loop here and add an 'area' field to the neurons in C_allspikes.
% I could also just treat C_allinfo like C_allspikes and flatten the cells
% across tetrodes into one array inside the main for-loop, then grab the
% area from C_allinfo. For now I am going to choose the first option.
C_allinfo = C_alldata(cellinfo_idx,:);

% Store all state data in one array. The order of the data
% stored in rows is that of stateNames.
C_allstates = C_alldata(states_idx,:);

% Linfield data
C_alllinf = C_alldata(linf_idx,:);

% Identify times at which velocity is above a threshold and fill a cell array
% full of structs that matches that of the other states like REM, NREM,
% etc. This is done by create_runstate.m.
pos_idx = find(contains(filetypes,'pos01'));
if isempty(pos_idx)
    error("pos01 data must be loaded to run this analysis.")
end
C_runstate = create_runstate(C_alldata(pos_idx,:));
% Do the same for times when vel < threshold.
C_stillstate = create_stillstate(C_alldata(pos_idx,:));

% Now add run and still state data to C_allstates 
C_allstates = [C_allstates; C_runstate; C_stillstate];
% for r = 1:size(C_alldata,2)
%     C_allstates{end+1,r} = C_runstate{1,r};
% end


% Add 'area' field to the spike data for later sorting. Also add place
% field data such as mean and std FR, sparsity, and spatial coverage for each trajectory.
% linfields must be loaded to add the placefield data.
for r = 1:size(C_allspikes,2)
    for e = 1:size(C_allspikes{1,r},2)
        for tet = 1:size(C_allspikes{1,r}{1,e},2)
            if ~isempty(C_allspikes{1,r}{1,e}{1,tet})
                for nrn = 1:size(C_allinfo{1,r}{1,e}{1,tet},2)
                    if ~isempty(C_allspikes{1,r}{1,e}{1,tet}{1,nrn})

                        % If the neuron exists in the spike file it should
                        % exist in the cellinfo file.
                        nrnCellinfo = C_allinfo{1,r}{1,e}{1,tet}{1,nrn};
                        

                        if isfield(nrnCellinfo,'area')
                            C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.area = nrnCellinfo.area; % Create field in spike struct.       
                        end
                        
                        % Linfields data is only available for behavioral
                        % epochs (even). Existence of neurons does not
                        % quite match that of C_allspikes.
                        if size(C_alllinf{1,r},2) >= e && ~isempty(C_alllinf{1,r}{1,e})
                            if ~isempty(C_alllinf{1,r}{1,e}{1,tet}) && size(C_alllinf{1,r}{1,e}{1,tet},2) >= nrn
                                if ~isempty(C_alllinf{1,r}{1,e}{1,tet}{1,nrn})
                                    trajsLinf = C_alllinf{1,r}{1,e}{1,tet}{1,nrn}; % Trajectories for a single neuron.
                                    for tr = 1:size(trajsLinf,2)
        
                                        trData = trajsLinf{tr};
                                        binRate = trData(:,7)./trData(:,6); % Firing rate normalized by occupancy
                                        FRmean = mean(binRate,1,'omitnan');
                                        FRstd = std(binRate,1,'omitnan');
        
                                        % calc sparsity
                                        totOcc = sum(trData(:,6),1); % Total time spent on track
                                        probOcc = trData(:,6)./totOcc; % Probability of occupying each spot on the trajectory
                                        sum1 = sum(probOcc.*binRate);
                                        sum2 = sum(probOcc.*(binRate.^2));
                                        sparsity = (sum1^2)/sum2;
        
                                        % calc spatial coverage
                                        numBins = size(trData,1); % number of spatial bins this trajectory is split into.
                                        peakRate = max(binRate);
                                        isAboveThr = binRate > 0.25*peakRate;
                                        coverage = sum(isAboveThr)/numBins; % Fraction of spatial bins with at least 25%
                                        % of the peak firing rate on that trajectory.

                                        
        
                                        % Assign to fields in struct
                                        C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_FRmeans(tr) = FRmean;
                                        C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_FRstds(tr) = FRstd;
                                        C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_sparsity(tr) = sparsity;
                                        C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_coverage(tr) = coverage;
                                        C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_FRpeaks(tr) = peakRate;
        
                                    end
                                    spStruct = C_allspikes{1,r}{1,e}{1,tet}{1,nrn};
                                    isPC = is_placecell(spStruct.traj_FRmeans, spStruct.traj_FRstds, spStruct.traj_FRpeaks,...
                                        spStruct.traj_sparsity, spStruct.traj_coverage);
                                    C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.isPlaceCell = isPC;

                                end
                            end
                        end

                    end
                end
            end
        end
    end
end




% ------------------------------------------------------------------

% Create firing rate matrices and place cell sum matrix

% Cell array to hold matrices for each brain region (in brainArea order).
% Matrices contain mean spike rates by epoch, not discriminating by state. 
% Rows are neuron identity, columns are epochs. This is for all rats.
FR_byEpoch = cell(1,length(brainAreas));

% Matrix to indicate how many trajectories qualify for place fields during that epoch, organized in the same format
% as FR_byEpoch. Each element of matrix contains an integer between 0 and 4.
PCsum_byEpoch = cell(1,length(brainAreas));


% Calculate and store mean firing rate across the entire epoch for each
% epoch, rat, and brain area. Same for isPC.
for r = 1:length(C_allspikes) 
    
    for a = 1:length(brainAreas) % I could eliminate this loop by using strfind(brainAreas) to index...

        % Creates a matrix to hold spike rate data for one rat. Dims (num nrns)
        % x (num epochs). The first epoch is used to find the number of nrns
        % because this does not change across epochs as stated above.
        ratRates = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));

        % Same as ratRates, but stores place cell identification info.
        ratPCsum = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));

        for e = 1:length(C_allspikes{1,r})
    
            nrnsAlltets = [C_allspikes{1,r}{1,e}{:}]; % Combining nrn data from all tets.
    
            for nrn = 1:length(nrnsAlltets)
                % isfield also returns false if the struct does not exist
                if isfield(nrnsAlltets{nrn},'area') && strcmp(nrnsAlltets{nrn}.area, brainAreas{a}) % If neuron is in wanted brain region

                    if isfield(nrnsAlltets{nrn},'data') 
                        % I checked and confirmed that calculating the mean firing rate 
                        % this way matches the .meanrate field value for each nrn (only off by 0.0001 for some nrns).
                        ratRates(nrn,e) = size(nrnsAlltets{nrn}.data,1)/(nrnsAlltets{nrn}.timerange(2)-nrnsAlltets{nrn}.timerange(1));
                    end
                
                    if isfield(nrnsAlltets{nrn},'isPlaceCell')
                        ratPCsum(nrn,e) = sum(nrnsAlltets{nrn}.isPlaceCell,2);
                    end

                end

            end
        end
        FR_byEpoch{a} = [FR_byEpoch{a}; ratRates]; 
        PCsum_byEpoch{a} = [PCsum_byEpoch{a}; ratPCsum];
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


% Create isPC_byEpoch, a similar matrix but filled with 1s and 0s denoting
% if the neuron is a place cell or not. Also, since only NaNs are in
% PCsum_byEpoch for rest epochs, a neuron in a rest epoch is considered a
% place cell if it was a place cell in the previous behavioral epoch. This
% is done here.
PCthr = 1;
for a = 1:length(brainAreas)
    isPC_byEpoch{a} = PCsum_byEpoch{a} >= PCthr;
    for col = 1:size(isPC_byEpoch{a},2)
        if ismember(col,restEpochs) && col > 1
            isPC_byEpoch{a}(:,col) = isPC_byEpoch{a}(:,col-1);
        end
    end
end


% Gets epochs that have data for all states. For most states this should be all odd
% epochs for the data I am sorting here, because all the states
% happen during the rest epochs except for rippletimes. I will detect non-empty
% epochs for generality. Note this assumes all rats have the same number of
% epochs as rat 1, which should be true when clip_17_epochs.m has been
% used on C_allstates.
nonemptyEpochs = zeros(size(C_allstates,1),length(C_allstates{1,1}));
for s2 = 1:size(C_allstates,1)
    nonemptyEpochs(s2,:) = ~cellfun(@isempty, C_allstates{s2,1});
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

    for s2 = 1:size(C_allstates,1)
    
        % Gets epochs that have data for this state. This should be identical across rats.
        nonempty = zeros(1,length(C_allstates{1,1}));
    
        for r = 1:size(C_allstates,2) 
            % Indicates which epochs contain data
            nonempty(1,:) = ~cellfun(@isempty, C_allstates{s2,r});
            % To hold mean firing rate data across all epochs for a single rat.
            ratMeans = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));
            % To hold state occurances for all epochs.
            occEpochs = cell(1,length(C_allstates{s2,r}));
    
            for e = 1:length(C_allstates{s2,r})
                
                if nonempty(1,e) == 1
                    % Now the real deal. For every non-empty epoch, each state type (rem, rippletime,
                    % sleep, sws, and waking) will loop through the times
                    % that that state occurs, counting the number of spikes
                    % that occur during the duration of that occurance. In the end
                    % these spike counts will be summed up and divided by the
                    % total duration of that state to get the mean rate.
                    
                    epochData = C_allstates{s2,r}{1,e}; % Data on a single rat, state, and epoch.
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
                    Dur_allStates(s2,r,e) = stateDur;
                    Dur_ratRef(s2,r,e) = r;
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

            FR_allStates{s2,a} = [FR_allStates{s2,a}; ratMeans];    
            Occ_allStates{s2,r} = occEpochs;
        end
    end
end




%% Plot histogram of velocities

% posData = C_alldata(find(contains(filetypes,'pos01')),:);
% allVelData = [];
% for r = 1:size(posData,2)
%     for e = 1:size(posData{1,r},2)
% 
%         allVelData = [allVelData; posData{1,r}{1,e}.data(:,5)];
% 
%     end
% end
% 
% edges = linspace(5, 55, sqrt(size(allVelData,1))); % There are very few velocity values above 55
% figure;
% histogram(allVelData, edges)
% title("Velocities Across all Rats and Epochs")
% ylabel("Count")
% xlabel("Velocity (cm/s)")


% %% I need to plot velocity with the times I have designated as "run" to check that 
% % what I did worked.
% 
% 
% posData = C_alldata(find(contains(filetypes,'pos01')),:); % Grabs position file data.
% runData = C_allstates(find(contains(stateNames,'run')),:);
% stillData = C_allstates(find(contains(stateNames,'still')),:);
% r = 1;
% e = 2;
% 
% posStruct = posData{1,r}{1,e};
% runStruct = runData{1,r}{1,e};
% stillStruct = stillData{1,r}{1,e};
% 
% vel = posStruct.data(:,5);
% vel_sm = sgolayfilt(vel, 2, 31); % Smoothing the data.
% 
% % timeRange = linspace(runStruct.timerange(1),runStruct.timerange(2),size(posStruct.data,1));
% 
% figure;
% hold on;
% title(sprintf("Velocity and Run State | Rat %d, Epoch %d", r, e))
% plot(posStruct.data(:,1), vel)
% plot(posStruct.data(:,1), vel_sm)
% 
% runLabels = {};
% for i = 1:length(runStruct.starttime)
%     x_vertices = [runStruct.starttime(i),runStruct.endtime(i),runStruct.endtime(i),runStruct.starttime(i)];
%     y_vertices = [min(vel_sm),min(vel_sm),max(vel_sm),max(vel_sm)];
%     patch(x_vertices, y_vertices, [0.6 0.9 0],'FaceAlpha', 0.3, 'EdgeColor','none')
% 
%     if i == 1
%         runLabels{i} = "run times";
%     else
%         runLabels{i} = "";
%     end
% end
% 
% stillLabels = {};
% for i = 1:length(stillStruct.starttime)
%     x_vertices = [stillStruct.starttime(i),stillStruct.endtime(i),stillStruct.endtime(i),stillStruct.starttime(i)];
%     y_vertices = [min(vel_sm),min(vel_sm),max(vel_sm),max(vel_sm)];
%     patch(x_vertices, y_vertices, [0 0.9 1],'FaceAlpha', 0.3, 'EdgeColor','none')
% 
%     if i == 1
%         stillLabels{i} = "still times";
%     else
%         stillLabels{i} = "";
%     end
% end
% 
% ylabel("Velocity (cm/s)")
% xlabel("Time (s)")
% legend({'vel','smooth vel', runLabels{:}, stillLabels{:}})





%% Begin FR scatter plot analysis
% NOTE this analysis treats cells between epochs as if they are
% different cells. I don't know if this is best, or if I should maintain
% cell identity across epochs and calculate an average FR.

% logicals to plot all individual plots and subplots.
plotindv = 0;
plotsubplots = 1;

stateColors = [[0 0.4470 0.7410]; [1 0.8 0]; [0.5 0.1 1]; [1 0 0]; [0.6 0.9 0]; [0 0.9 1]];
    
numg = 1; % Number of groups to evenly split the data into.


while(1) % Prompt user to determine whether or not to continue with plotting.
    usri = input(sprintf('%d figures will be created. Do you wish to continue? y/n :',...
        plotindv*(size(FR_allStates,1)^2)*length(brainAreas) + plotsubplots*size(FR_allStates,1)*length(brainAreas)),'s');
    if strcmpi(usri,'y')
        disp("Plotting figures...")
        break;
    else
        disp("Execution stopped.")
        return
    end
end


%I will loop through all the states twice. This way I plot all
%combinations of states vs. eachother.
for s1 = 1:size(FR_allStates,1)
    
    if plotsubplots

        for a = 1:length(brainAreas)
            figure;
            sgtitle(sprintf("Neuron FRs each Epoch| %s vs. all States | %s",stateNames{s1},brainAreas{a}))
            for s2 = 1:size(FR_allStates,1)
                
                % Form vectors of firing rates for the two states being plotted against
                % each other, making sure to perserve the cell identity in each vector.
                % Each vector should hold info from all epochs.                
                s1Data = FR_allStates{s1,a};
                s2Data = FR_allStates{s2,a}; % Data of state to be plotted against.

                % For run and still states I only want the behavioral
                % epochs. All other states I want rest epochs.
                if any(contains(['run','still'],stateNames{s1})) && any(contains(['run','still'],stateNames{s2}))
                    s1Data = s1Data(:,behEpochs);
                    s2Data = s2Data(:,behEpochs);
                elseif any(contains(['run','still'],stateNames{s1}))
                    s1Data = s1Data(:,behEpochs);
                    % Temporay fix. Omits first sleep epoch to keep number of epochs the same between behavior and rest.
                    s2Data = s2Data(:,restEpochs(2:end)); 
                elseif any(contains(['run','still'],stateNames{s2}))
                    s1Data = s1Data(:,restEpochs(2:end));
                    s2Data = s2Data(:,behEpochs);
                else
                    s1Data = s1Data(:,restEpochs);
                    s2Data = s2Data(:,restEpochs);
                end
            
                
                % NOTE I need to leave the NaNs in these matrices to preserve the
                % neuron IDs. Some neurons have NaN in one matrix and not the
                % other, so removing NaNs messes up the matching of x and y values
                % between the two lists. 
                
                subplot(ceil(size(FR_allStates,1)/2),2,s2);
                hold on;
                groups = randi([1,numg],[length(s1Data(:)),1]);
                clrs = [];
                if numg > 1
                    for g = 1:numg-1
                        clrs = [clrs; [0.7 0.7 0.7]]; % Assign all but last group to color grey
                    end
                end
                clrs = [clrs; stateColors(s2,:)]; 
        
                gscatter(gca, s1Data(:),s2Data(:),groups, clrs, '.', 1, doleg='off');
                plot(min(s2Data(:)):max(s2Data(:)), min(s2Data(:)):max(s2Data(:)),'k')
                xlabel(sprintf("%s Mean FR (Hz)",stateNames{s1}))
                ylabel(sprintf("%s Mean FR (Hz)", stateNames{s2}))
                set(gca, 'XScale', 'log', 'YScale', 'log');
            end
        end
    end
    
    
    if plotindv % Individual FR plots
    
        for a = 1:length(brainAreas)
            for s2 = 1:size(FR_allStates,1)
                
                % Form vectors of firing rates for the two states being plotted against
                % each other, making sure to perserve the cell identity in each vector.
                % Each vector should hold info from all epochs.                
                s1Data = FR_allStates{s1,a};
                s2Data = FR_allStates{s2,a}; % Data of state to be plotted against.
            
                s1Data = s1Data(:,restEpochs);
                s2Data = s2Data(:,restEpochs);
                % NOTE I need to leave the NaNs in these matrices to preserve the
                % neuron IDs. Some neurons have NaN in one matrix and not the
                % other, so removing NaNs messes up the matching of x and y values
                % between the two lists. 
                
                figure;
                hold on;
                groups = randi([1,numg],[length(s1Data(:)),1]);
                clrs = [];
                if numg > 1
                    for g = 1:numg-1
                        clrs = [clrs; [0.7 0.7 0.7]]; % Assign all but last group to color grey
                    end
                end
                clrs = [clrs; stateColors(s2,:)]; 
        
                gscatter(gca, s1Data(:),s2Data(:),groups, clrs, '.', 1, doleg='off');
                plot(min(s2Data(:)):max(s2Data(:)), min(s2Data(:)):max(s2Data(:)),'k')
                xlabel(sprintf("%s Mean FR (Hz)",stateNames{s1}))
                ylabel(sprintf("%s Mean FR (Hz)", stateNames{s2}))
                title(sprintf("%s vs. %s | %s",stateNames{s1},stateNames{s2},brainAreas{a}))
                set(gca, 'XScale', 'log', 'YScale', 'log');
            end
        end
    
    end

end



%% Identify place cells and interneurons and label them on the scatter plots

% How blake identifies place cells:
% He did a lot of research and found that there is no concensus on how to
% define place cells. How he does it:
% Cells must satisfy threshold values for these measures:
% Mean firing rate - below 10 Hz or so (double check).
% Peak firing rate - Must be x number of std above mean FR and above 3 Hz.
% Number of spikes - More than 100 spikes.
% Sparsity - From my observations, less than 0.5 seems good
% Spatial specificity
% Spatial coherence
% Spatial coverage - As defined in Jadhav, Frank 2016: "the fraction of
% area above 25% of the peak firing rate across all four trajectories for
% neurons with a peak spatial firing rate > 3Hz."

% NOTE: When I combine all trajectories into one, the sparsity and spatial
% coverage metrics can break down and provide very low values for cells
% that are clearly not place cells.



% % Combining all trajectories into one and running analyses on it. Sparsity
% % and spatial coverage seem to break down.
% for r = 1:size(C_alllinf,2)
%     for e = 1:size(C_alllinf{1,r},2)
% 
%         if ~isempty(C_alllinf{1,r}{1,e}) % Goes into behavioral epochs only.
% 
%             nrnsAllTets = [C_alllinf{1,r}{1,e}{:}]; % Gets all neurons for one rat and one epoch.
% 
%             for nrn = 1:size(nrnsAllTets,2)
% 
%                 if ~isempty(nrnsAllTets{1,nrn})
%                     figure;
%                     title(sprintf("Epoch %d Nrn %d, all trajs | spar=%.2f, cov=%.2f",e,nrn,sparsity,coverage))
% 
%                     catTrData = vertcat(nrnsAllTets{1,nrn}{:}); % Concatenated data across all 4 trajectories.
%                     FRmean = mean(catTrData(:,5),1,'omitnan');
%                     FRstd = std(catTrData(:,5),1,'omitnan');
% 
%                     % calc sparsity
%                     totOcc = sum(catTrData(:,6),1); % Total time spent on track
%                     probOcc = catTrData(:,6)./totOcc; % Probability of occupying each spot on the trajectory
%                     binRate = catTrData(:,7)./catTrData(:,6); % Firing rate normalized by occupancy
%                     sum1 = sum(probOcc.*binRate);
%                     sum2 = sum(probOcc.*(binRate.^2));
%                     sparsity = (sum1^2)/sum2;
% 
%                     % calc spatial coverage
%                     numBins = size(catTrData,1); % number of spatial bins this trajectory is split into.
%                     peakRate = max(binRate);
%                     isAboveThr = binRate > 0.25*peakRate;
%                     coverage = sum(isAboveThr)/numBins; % Fraction of spatial bins with at least 25%
%                     % of the peak firing rate on that trajectory.
% 
%                     xpos = 1:size(catTrData,1);
%                     xpos = (xpos-1).*2; % Transform into cm
% 
%                     hold on;
%                     plot(xpos, binRate, 'k')
%                     plot(xpos, ones(size(xpos))*FRmean, "--r")
%                     plot(xpos, ones(size(xpos))*(FRmean+FRstd), "--b")
%                     plot(xpos, ones(size(xpos))*(FRmean+2*FRstd), "--b")
%                     plot(xpos, ones(size(xpos))*(FRmean+3*FRstd), "--b")
%                     ylabel("Occ Norm FR")
%                     xlabel("Lin Pos All Trajs (cm)")
% 
%                     pause
%                     close
%                 end
% 
% 
%             end
% 
%         end
% 
%     end
% end


% % Analyzing all trajectories seperately.
% for r = 1:size(C_alllinf,2)
%     for e = 1:size(C_alllinf{1,r},2)
% 
%         if ~isempty(C_alllinf{1,r}{1,e}) % Goes into behavioral epochs only.
% 
%             nrnsAllTets = [C_alllinf{1,r}{1,e}{:}]; % Gets all neurons for one rat and one epoch.
% 
%             for nrn = 1:size(nrnsAllTets,2)
% 
%                 if ~isempty(nrnsAllTets{1,nrn})
%                     figure;
%                     sgtitle(sprintf("Epoch %d Nrn %d",e,nrn))
%                     for tr = 1:4 % 4 different trajectories
%                         trData = nrnsAllTets{1,nrn}{1,tr};
%                         FRmean = mean(trData(:,5),1,'omitnan');
%                         FRstd = std(trData(:,5),1,'omitnan');
% 
%                         % calc sparsity
%                         totOcc = sum(trData(:,6),1); % Total time spent on track
%                         probOcc = trData(:,6)./totOcc; % Probability of occupying each spot on the trajectory
%                         binRate = trData(:,7)./trData(:,6); % Firing rate normalized by occupancy
%                         sum1 = sum(probOcc.*binRate);
%                         sum2 = sum(probOcc.*(binRate.^2));
%                         sparsity = (sum1^2)/sum2;
% 
%                         % calc spatial coverage
%                         numBins = size(trData,1); % number of spatial bins this trajectory is split into.
%                         peakRate = max(binRate);
%                         isAboveThr = binRate > 0.25*peakRate;
%                         coverage = sum(isAboveThr)/numBins; % Fraction of spatial bins with at least 25%
%                         % of the peak firing rate on that trajectory.
% 
%                         ax(tr) = subplot(2,2,tr);
%                         hold on;
%                         plot(trData(:,1), binRate, 'k')
%                         plot(trData(:,1), ones(size(trData(:,1)))*FRmean, "--r")
%                         plot(trData(:,1), ones(size(trData(:,1)))*(FRmean+FRstd), "--b")
%                         plot(trData(:,1), ones(size(trData(:,1)))*(FRmean+2*FRstd), "--b")
%                         plot(trData(:,1), ones(size(trData(:,1)))*(FRmean+3*FRstd), "--b")
%                         ylabel("Occ Norm FR")
%                         xlabel("Lin Pos (cm)")
%                         title(sprintf("Traj %d | spar=%.2f, cov=%.2f",tr,sparsity,coverage))
%                         linkaxes(ax,'y')
%                     end
%                     pause
%                     close
%                 end
% 
% 
%             end
% 
%         end
% 
%     end
% end


% % Compare occupancy and spike count.
% for r = 1:size(C_alllinf,2)
%     for e = 1:size(C_alllinf{1,r},2)
% 
%         if ~isempty(C_alllinf{1,r}{1,e}) % Goes into behavioral epochs only.
% 
%             nrnsAllTets = [C_alllinf{1,r}{1,e}{:}]; % Gets all neurons for one rat and one epoch.
% 
%             for nrn = 1:size(nrnsAllTets,2)
% 
%                 if ~isempty(nrnsAllTets{1,nrn})
% 
% 
%                     for tr = 1:4 % 4 different trajectories
%                         figure;
%                         sgtitle(sprintf("Epoch %d Nrn %d Traj %d",e,nrn,tr))
%                         trData = nrnsAllTets{1,nrn}{1,tr};
%                         FRmean = mean(trData(:,5),1,'omitnan');
%                         FRstd = std(trData(:,5),1,'omitnan');
% 
%                         subplot(2,2,1);
%                         hold on;
%                         plot(trData(:,1), trData(:,5), 'k')
%                         plot(trData(:,1), ones(size(trData(:,1)))*FRmean, "--r")
%                         plot(trData(:,1), ones(size(trData(:,1)))*(FRmean+FRstd), "--b")
%                         plot(trData(:,1), ones(size(trData(:,1)))*(FRmean+2*FRstd), "--b")
%                         plot(trData(:,1), ones(size(trData(:,1)))*(FRmean+3*FRstd), "--b")
%                         ylabel("Occ Norm FR")
%                         xlabel("Lin Pos (cm)")
%                         title("Occ Norm FR")
% 
%                         subplot(2,2,2)
%                         hold on;
%                         plot(trData(:,1), trData(:,6), 'b')
%                         xlabel("Lin Pos (cm)")
%                         title("Smooth Occ")
% 
%                         subplot(2,2,3)
%                         hold on;
%                         plot(trData(:,1), trData(:,7), 'r')
%                         xlabel("Lin Pos (cm)")
%                         title("Smooth Spike Count")
% 
%                         subplot(2,2,4)
%                         hold on;
%                         plot(trData(:,1), trData(:,7)./trData(:,6))
%                         xlabel("Lin Pos (cm)")
%                         title("Occ/SC")
% 
%                         pause
%                         close
% 
%                     end
% 
%                 end
% 
% 
%             end
% 
%         end
% 
%     end
% end


% % Plot distribution of cells as the number of trajectories they qualify as
% % place cells on.
% pcSums = [];
% for r = 1:size(C_alllinf,2)
%     for e=2:2:size(C_alllinf{1,r},2)
%         for tet=1:size(C_alllinf{1,r}{1,e},2)
%             if ~isempty(C_alllinf{1,r}{1,e}{1,tet})
%                 for nrn = 1:size(C_alllinf{1,r}{1,e}{1,tet},2)
%                     if ~isempty(C_alllinf{1,r}{1,e}{1,tet}{1,nrn}) && ~isempty(C_allspikes{1,r}{1,e}{1,tet}{1,nrn})
% 
%                         spStruct = C_allspikes{1,r}{1,e}{1,tet}{1,nrn};
%                         pcSums(end+1) = sum(spStruct.isPlaceCell);
% 
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% figure;
% histogram(pcSums)
% title("Number of Place Cell-Qualifying Trajectories")
% ylabel("Count")
% xlabel("Number of Trajectories")


% % Plot how the number of place cell-qualifying trajectories changes over
% % the epochs for each neuron.
% for a = 1:size(PCsum_byEpoch,2) % Loop both brain areas
% 
    % PCsumMean = mean(PCsum_byEpoch{a}(:,behEpochs),1,'omitnan');
    % PCsumSTD = std(PCsum_byEpoch{a}(:,behEpochs),[],1,'omitnan');
    % PCsumSEM = PCsumSTD./sqrt(sum(~isnan(PCsum_byEpoch{a}(:,behEpochs)),1));
    % 
    % figure;
    % title(sprintf("Number of Place Cell-Qualifying Trajectories | %s",brainAreas{a}))
    % hold on;
    % ylabel("Number of Trajectories")
    % xlabel("Epoch")
    % shadedErrorBar(behEpochs, PCsumMean, PCsumSEM)
% 
% end

% % Plot this all together
% for a = 1:size(PCsum_byEpoch,2)
%     figure;
%     hold on;
%     title(sprintf("Fraction of Place Cells by Varying T | %s",brainAreas{a}));
%     ylabel("Fraction of Cells")
%     xlabel("Epoch")
%     legend()
% 
%     for PCthr = 0:4
%         isPCmat = PCsum_byEpoch{a}(:,behEpochs) >= PCthr;
%         plot(behEpochs, sum(isPCmat,1)./max(sum(isPCmat,1)), DisplayName = num2str(PCthr))
% 
%     end
% end



% Now actually plot new scatter plots

% stateColors = [[0 0.4470 0.7410]; [1 0.8 0]; [0.5 0.1 1]; [1 0 0]; [0.6 0.9 0]; [0 0.9 1]];
stateColors = [ [0 0.4470 0.7410]; [1 0 0]; [0.5 0.1 1]; [0.6 0.9 0]; [0 0.9 1]];


% while(1) % Prompt user to determine whether or not to continue with plotting.
%     usri = input(sprintf('%d figures will be created. Do you wish to continue? y/n :',...
%         plotindv*(size(FR_allStates,1)^2)*length(brainAreas) + plotsubplots*size(FR_allStates,1)*length(brainAreas)),'s');
%     if strcmpi(usri,'y')
%         disp("Plotting figures...")
%         break;
%     else
%         disp("Execution stopped.")
%         return
%     end
% end

% A note on my scatter plots: The matrices containing the data have epochs
% across the columns, and rows are concatenated neuron data for all rats.
% However, when I plot this data, I flatten it all into one column vector.
% This means that the scatter plot plots one point for every neuron from
% every epoch from every rat. So the same neuron will appear multiple times
% in the plot, once for every epoch. These scatter plots compare across
% states, so I want to compare the same neuron's firing rate across different
% states. Sometimes I compare behavioral and rest epochs, so then I am not
% comparing neurons from the same epochs, but rather from adjacent behavior and rest epochs.
for s1 = 1:size(FR_allStates,1)

    for a = 1:length(brainAreas)
        
        
        figure;
        sgtitle(sprintf("Neuron FRs each Epoch (black=PC, brown=int) | %s vs. all States | %s",stateNames{s1},brainAreas{a}))
        for s2 = 1:size(FR_allStates,1)
            
            % Form vectors of firing rates for the two states being plotted against
            % each other, making sure to perserve the cell identity in each vector.
            % Each vector should hold info from all epochs.                
            s1Data = FR_allStates{s1,a};
            s2Data = FR_allStates{s2,a}; % Data of state to be plotted against.



            % For run and still states I only want the behavioral
            % epochs. All other states I want rest epochs.
            if any(contains(['run','still'],stateNames{s1})) && any(contains(['run','still'],stateNames{s2}))
                s1Data = s1Data(:,behEpochs);
                s2Data = s2Data(:,behEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s1eFRmat = FR_byEpoch{1,a}(:,behEpochs);
                s2eFRmat = FR_byEpoch{1,a}(:,behEpochs);
            elseif any(contains(['run','still'],stateNames{s1}))
                s1Data = s1Data(:,behEpochs);
                % Temporay fix. Omits first sleep epoch to keep number
                % of epochs the same between behavior and rest. This
                % means that each data point is comparing the behavior
                % firing rate to the sleep firing rate in the next rest
                % epoch.
                s2Data = s2Data(:,restEpochs(2:end));
                s1isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,restEpochs(2:end));
                s1eFRmat = FR_byEpoch{1,a}(:,behEpochs);
                s2eFRmat = FR_byEpoch{1,a}(:,restEpochs(2:end));
            elseif any(contains(['run','still'],stateNames{s2}))
                s1Data = s1Data(:,restEpochs(2:end));
                s2Data = s2Data(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,restEpochs(2:end));
                s2eFRmat = FR_byEpoch{1,a}(:,behEpochs);
                s1eFRmat = FR_byEpoch{1,a}(:,restEpochs(2:end));
            else
                s1Data = s1Data(:,restEpochs);
                s2Data = s2Data(:,restEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,restEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,restEpochs);
                s1eFRmat = FR_byEpoch{1,a}(:,restEpochs);
                s2eFRmat = FR_byEpoch{1,a}(:,restEpochs);
            end
            
            
            % I want to separate place cells and non-PCs. But, to
            % maintain neuronal identity in the matrices, I can't just
            % remove one category or the other from sxData because
            % there can be different numbers of place cells between two
            % states (actually, because whether or not a cell is
            % considered a PC during a rest epoch is taken from the
            % previous behavioral epoch, and because in the scatter
            % plots behavioral epochs are compared to the subsequent
            % rest epoch, it would probably work to just remove the
            % unwanted neurons from the matrices before plotting).
            % However, to be safe, and in case this ever changes, I
            % want to replace unwanted neurons in sxData with NaN to be
            % sure that no misalignment happens between s1Data and
            % s2Data.
            % Create matrices with only PC data.
            s1DataPCs = s1Data;
            s1DataPCs(~s1isPCmat) = NaN;
            s2DataPCs = s2Data;
            s2DataPCs(~s2isPCmat) = NaN;
            % Create matrices with only non-PC data.
            s1DataNonPCs = s1Data;
            s1DataNonPCs(s1isPCmat) = NaN;
            s2DataNonPCs = s2Data;
            s2DataNonPCs(s2isPCmat) = NaN;

            % Label interneurons as those with a mean FR across the entire
            % epoch above a threshold intThr Hz.
            intThr = 7;
            s1isInt = s1eFRmat > intThr;
            s2isInt = s2eFRmat > intThr;
            % Also only label interneurons from the non-PC pool.
            s1DataInts = s1DataNonPCs;
            s2DataInts = s2DataNonPCs;
            s1DataInts(~s1isInt) = NaN;
            s2DataInts(~s2isInt) = NaN;


            % Note that plot only includes data points that have no NaN
            % values in x or y, so if a neuron/epoch is only active in one
            % state and not the other, this neuron/epoch is not plotted.
            subplot(ceil(size(FR_allStates,1)/2),2,s2);
            hold on;
            plot(s1DataNonPCs(:),s2DataNonPCs(:), '.', Color=stateColors(s2,:), MarkerSize=2);
            plot(s1DataPCs(:),s2DataPCs(:), '.', Color='k', MarkerSize=2);
            plot(s1DataInts(:),s2DataInts(:), '.', Color=[1 0.4134 0.1788], MarkerSize=2)
            plot(min(s2Data(:)):max(s2Data(:)), min(s2Data(:)):max(s2Data(:)),'k')
            xlabel(sprintf("%s Mean FR (Hz)",stateNames{s1}))
            ylabel(sprintf("%s Mean FR (Hz)", stateNames{s2}))
            set(gca, 'XScale', 'log', 'YScale', 'log');
        end
    end

end





%% Plot place cells color coded by spatial coverage

% I want to try and bring the information from the neurons in C_allspikes
% in a different way than I did for the isPC information above. This time,
% I want to loop through all neurons of all rats, effectively going row by row
% in a FR_allStates submatrix. This way I can keep track of which neuron has
% which properties and easily add them to the plot without having to make a
% separate matrix all the way back when FR_byEpoch is created...
% But... I'm not sure how to do this. In the end I need a matrix or array
% that has the coverage for each neuron in the same order as FR_byEpoch.
% So, until I can think of a better way, I am going to create the matrix in
% the same way as FR_byEpoch.

% Because there are four measures of sparsity (one for each traj) for each
% neuron, I am going to use isPlaceCell to tell me which trajs a neuron has
% a place field on, and then average the spatial coverage from all those
% trajs to get a single measure of coverage for each neuron.
cov_byEpoch = cell(1,length(brainAreas));

% Do the exact same thing, but for peak firing rate on all place field
% trajs.
peak_byEpoch = cell(1,length(brainAreas));

for r = 1:size(C_allspikes,2)

    for a = 1:length(brainAreas)

        ratCovs = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));
        ratPeaks = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));

        for e = 1:size(C_allspikes{1,r},2)
            
            nrnsAlltets = [C_allspikes{1,r}{1,e}{:}]; % Combining nrn data from all tets.
    
            for nrn = 1:length(nrnsAlltets)
                % isfield also returns false if the struct does not exist
                if isfield(nrnsAlltets{nrn},'isPlaceCell') && strcmp(nrnsAlltets{nrn}.area, brainAreas{a}) % If neuron is in wanted brain region
                
                    % Multiply coverage by isPlaceCell to only get spatial
                    % coverage of trajectories that have place fields.
                    PCcov = nonzeros(nrnsAlltets{nrn}.isPlaceCell.*nrnsAlltets{nrn}.traj_coverage);
                    ratCovs(nrn,e) = mean(PCcov); % The mean coverage from all place field trajs.

                    % Peaks of trajs with place fields
                    PCpeak = nonzeros(nrnsAlltets{nrn}.isPlaceCell.*nrnsAlltets{nrn}.traj_FRpeaks);
                    ratPeaks(nrn,e) = mean(PCpeak);

                end

            end
        end

        cov_byEpoch{a} = [cov_byEpoch{a}; ratCovs];
        peak_byEpoch{a} = [peak_byEpoch{a}; ratPeaks];
    end
end

% Fill in rest epoch with previous beh epoch data.
for a = 1:length(brainAreas)
    for col = 1:size(cov_byEpoch{a},2)
        if ismember(col,restEpochs) && col > 1
            cov_byEpoch{a}(:,col) = cov_byEpoch{a}(:,col-1);
            peak_byEpoch{a}(:,col) = peak_byEpoch{a}(:,col-1);
        end
    end
end




% Scatter plots for spatial coverage

for s1 = 1:size(FR_allStates,1)

    for a = 1:length(brainAreas)


        figure;
        sgtitle(sprintf("Place Cell FRs and Coverage for each Epoch  | %s vs. all States | %s",stateNames{s1},brainAreas{a}))
        for s2 = 1:size(FR_allStates,1)

            % Form vectors of firing rates for the two states being plotted against
            % each other, making sure to perserve the cell identity in each vector.
            % Each vector should hold info from all epochs.                
            s1Data = FR_allStates{s1,a};
            s2Data = FR_allStates{s2,a}; % Data of state to be plotted against.

            % For run and still states I only want the behavioral
            % epochs. All other states I want rest epochs.
            if any(contains(['run','still'],stateNames{s1})) && any(contains(['run','still'],stateNames{s2}))
                s1Data = s1Data(:,behEpochs);
                s2Data = s2Data(:,behEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s1covmat = cov_byEpoch{1,a}(:,behEpochs);
                s2covmat = cov_byEpoch{1,a}(:,behEpochs);
            elseif any(contains(['run','still'],stateNames{s1}))
                s1Data = s1Data(:,behEpochs);
                s2Data = s2Data(:,restEpochs(2:end));
                s1isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,restEpochs(2:end));
                s1covmat = cov_byEpoch{1,a}(:,behEpochs);
                s2covmat = cov_byEpoch{1,a}(:,restEpochs(2:end));
            elseif any(contains(['run','still'],stateNames{s2}))
                s1Data = s1Data(:,restEpochs(2:end));
                s2Data = s2Data(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,restEpochs(2:end));
                s2covmat = cov_byEpoch{1,a}(:,behEpochs);
                s1covmat = cov_byEpoch{1,a}(:,restEpochs(2:end));
            else
                s1Data = s1Data(:,restEpochs);
                s2Data = s2Data(:,restEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,restEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,restEpochs);
                s1covmat = cov_byEpoch{1,a}(:,restEpochs);
                s2covmat = cov_byEpoch{1,a}(:,restEpochs);
            end

            % Create matrices with only PC data.
            s1DataPCs = s1Data;
            s1DataPCs(~s1isPCmat) = NaN;
            s2DataPCs = s2Data;
            s2DataPCs(~s2isPCmat) = NaN;

            % The fact that a neuron has a value in cov_byEpoch means
            % that it has at least one traj with a place field, but we
            % want to be able to adjust the traj threshold PCthr, so we
            % will filter the coverage data by the isPC data.
            s1covmat(~s1isPCmat) = NaN;
            s2covmat(~s2isPCmat) = NaN;

            % The spatial coverage valuyes should be the same between
            % the two states, but lets average just in case.
            catMat = cat(3,s1covmat,s2covmat);
            colorCovMat = mean(catMat,3);


            subplot(ceil(size(FR_allStates,1)/2),2,s2);
            hold on;
            scatter(s1DataPCs(:),s2DataPCs(:), 30, colorCovMat(:),'.');
            c = colorbar;
            clim([0,0.5]);
            colormap jet
            plot(min(s2DataPCs(:)):max(s2DataPCs(:)), min(s2DataPCs(:)):max(s2DataPCs(:)),'k')
            xlabel(sprintf("%s Mean FR (Hz)",stateNames{s1}))
            ylabel(sprintf("%s Mean FR (Hz)", stateNames{s2}))
            set(gca, 'XScale', 'log', 'YScale', 'log');
        end
    end
end




%% Plot peak firing rate as color code across FR scatter plot of only place cells

% Scatter plots for peak FR
for s1 = 1:size(FR_allStates,1)

    for a = 1:length(brainAreas)


        figure;
        sgtitle(sprintf("Place Cell Peak FR for each Epoch  | %s vs. all States | %s",stateNames{s1},brainAreas{a}))
        for s2 = 1:size(FR_allStates,1)

            % Form vectors of firing rates for the two states being plotted against
            % each other, making sure to perserve the cell identity in each vector.
            % Each vector should hold info from all epochs.                
            s1Data = FR_allStates{s1,a};
            s2Data = FR_allStates{s2,a}; % Data of state to be plotted against.

            % For run and still states I only want the behavioral
            % epochs. All other states I want rest epochs.
            if any(contains(['run','still'],stateNames{s1})) && any(contains(['run','still'],stateNames{s2}))
                s1Data = s1Data(:,behEpochs);
                s2Data = s2Data(:,behEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s1peakmat = peak_byEpoch{1,a}(:,behEpochs);
                s2peakmat = peak_byEpoch{1,a}(:,behEpochs);
            elseif any(contains(['run','still'],stateNames{s1}))
                s1Data = s1Data(:,behEpochs);
                s2Data = s2Data(:,restEpochs(2:end));
                s1isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,restEpochs(2:end));
                s1peakmat = peak_byEpoch{1,a}(:,behEpochs);
                s2peakmat = peak_byEpoch{1,a}(:,restEpochs(2:end));
            elseif any(contains(['run','still'],stateNames{s2}))
                s1Data = s1Data(:,restEpochs(2:end));
                s2Data = s2Data(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,restEpochs(2:end));
                s2peakmat = peak_byEpoch{1,a}(:,behEpochs);
                s1peakmat = peak_byEpoch{1,a}(:,restEpochs(2:end));
            else
                s1Data = s1Data(:,restEpochs);
                s2Data = s2Data(:,restEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,restEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,restEpochs);
                s1peakmat = peak_byEpoch{1,a}(:,restEpochs);
                s2peakmat = peak_byEpoch{1,a}(:,restEpochs);
            end

            % Create matrices with only PC data.
            s1DataPCs = s1Data;
            s1DataPCs(~s1isPCmat) = NaN;
            s2DataPCs = s2Data;
            s2DataPCs(~s2isPCmat) = NaN;

            s1peakmat(~s1isPCmat) = NaN;
            s2peakmat(~s2isPCmat) = NaN;

            catMat = cat(3,s1peakmat,s2peakmat);
            colorPeakMat = mean(catMat,3);


            subplot(ceil(size(FR_allStates,1)/2),2,s2);
            hold on;
            scatter(s1DataPCs(:),s2DataPCs(:), 20, colorPeakMat(:),'.');
            c = colorbar;
            colormap jet
            clim([0,30])
            plot(min(s2DataPCs(:)):max(s2DataPCs(:)), min(s2DataPCs(:)):max(s2DataPCs(:)),'k')
            xlabel(sprintf("%s Mean FR (Hz)",stateNames{s1}))
            ylabel(sprintf("%s Mean FR (Hz)", stateNames{s2}))
            set(gca, 'XScale', 'log', 'YScale', 'log');
        end
    end
end



%% Plot fold change in FR with spatial coverage and peak FR

% Of course this can only be done for place cells.

cov_byEpoch = cell(1,length(brainAreas));
peak_byEpoch = cell(1,length(brainAreas));

for r = 1:size(C_allspikes,2)

    for a = 1:length(brainAreas)

        ratCovs = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));
        ratPeaks = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));

        for e = 1:size(C_allspikes{1,r},2)
            
            nrnsAlltets = [C_allspikes{1,r}{1,e}{:}]; % Combining nrn data from all tets.
    
            for nrn = 1:length(nrnsAlltets)
                % isfield also returns false if the struct does not exist
                if isfield(nrnsAlltets{nrn},'isPlaceCell') && strcmp(nrnsAlltets{nrn}.area, brainAreas{a}) % If neuron is in wanted brain region
                
                    % Multiply coverage by isPlaceCell to only get spatial
                    % coverage of trajectories that have place fields.
                    PCcov = nonzeros(nrnsAlltets{nrn}.isPlaceCell.*nrnsAlltets{nrn}.traj_coverage);
                    ratCovs(nrn,e) = mean(PCcov); % The mean coverage from all place field trajs.

                    % Peaks of trajs with place fields
                    PCpeak = nonzeros(nrnsAlltets{nrn}.isPlaceCell.*nrnsAlltets{nrn}.traj_FRpeaks);
                    ratPeaks(nrn,e) = mean(PCpeak);

                end

            end
        end

        cov_byEpoch{a} = [cov_byEpoch{a}; ratCovs];
        peak_byEpoch{a} = [peak_byEpoch{a}; ratPeaks];
    end
end

% Fill in rest epoch with previous beh epoch data.
for a = 1:length(brainAreas)
    for col = 1:size(cov_byEpoch{a},2)
        if ismember(col,restEpochs) && col > 1
            cov_byEpoch{a}(:,col) = cov_byEpoch{a}(:,col-1);
            peak_byEpoch{a}(:,col) = peak_byEpoch{a}(:,col-1);
        end
    end
end



% Scatter plots for spatial coverage

% stateColors = [[0 0.4470 0.7410]; [1 0.8 0]; [0.5 0.1 1]; [1 0 0]; [0.6 0.9 0]; [0 0.9 1]];
stateColors = [ [0 0.4470 0.7410]; [1 0 0]; [0.5 0.1 1]; [0.6 0.9 0]; [0 0.9 1]];

for s1 = 1:size(FR_allStates,1)

    for a = 1:length(brainAreas)


        figure;
        sgtitle(sprintf("Place Cell Coverage vs Fold Change in FR| %s vs. all States | %s",stateNames{s1},brainAreas{a}))
        for s2 = 1:size(FR_allStates,1)

            % Form vectors of firing rates for the two states being plotted against
            % each other, making sure to perserve the cell identity in each vector.
            % Each vector should hold info from all epochs.                
            s1Data = FR_allStates{s1,a};
            s2Data = FR_allStates{s2,a}; % Data of state to be plotted against.

            % For run and still states I only want the behavioral
            % epochs. All other states I want rest epochs.
            if any(contains(['run','still'],stateNames{s1})) && any(contains(['run','still'],stateNames{s2}))
                s1Data = s1Data(:,behEpochs);
                s2Data = s2Data(:,behEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s1covmat = cov_byEpoch{1,a}(:,behEpochs);
                s2covmat = cov_byEpoch{1,a}(:,behEpochs);
            elseif any(contains(['run','still'],stateNames{s1}))
                s1Data = s1Data(:,behEpochs);
                s2Data = s2Data(:,restEpochs(2:end));
                s1isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,restEpochs(2:end));
                s1covmat = cov_byEpoch{1,a}(:,behEpochs);
                s2covmat = cov_byEpoch{1,a}(:,restEpochs(2:end));
            elseif any(contains(['run','still'],stateNames{s2}))
                s1Data = s1Data(:,restEpochs(2:end));
                s2Data = s2Data(:,behEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,behEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,restEpochs(2:end));
                s2covmat = cov_byEpoch{1,a}(:,behEpochs);
                s1covmat = cov_byEpoch{1,a}(:,restEpochs(2:end));
            else
                s1Data = s1Data(:,restEpochs);
                s2Data = s2Data(:,restEpochs);
                s1isPCmat = isPC_byEpoch{1,a}(:,restEpochs);
                s2isPCmat = isPC_byEpoch{1,a}(:,restEpochs);
                s1covmat = cov_byEpoch{1,a}(:,restEpochs);
                s2covmat = cov_byEpoch{1,a}(:,restEpochs);
            end

            % Create matrices with only PC data.
            s1DataPCs = s1Data;
            s1DataPCs(~s1isPCmat) = NaN;
            s2DataPCs = s2Data;
            s2DataPCs(~s2isPCmat) = NaN;

            % The fact that a neuron has a value in cov_byEpoch means
            % that it has at least one traj with a place field, but we
            % want to be able to adjust the traj threshold PCthr, so we
            % will filter the coverage data by the isPC data.
            s1covmat(~s1isPCmat) = NaN;
            s2covmat(~s2isPCmat) = NaN;

            % The spatial coverage values should be the same between
            % the two states, but lets average just in case.
            catMat = cat(3,s1covmat,s2covmat);
            colorCovMat = mean(catMat,3);
            


            s2s1fold = (s2DataPCs-s1DataPCs)./s1DataPCs; % Fold change in state 2 FR with respect to state 1 FR.


            subplot(ceil(size(FR_allStates,1)/2),2,s2);
            hold on;
            % plot(s2s1fold(:), colorCovMat(:), '.')

            scatter(s2s1fold(:),colorCovMat(:), 20, stateColors(s2,:),'.');
            % c = colorbar;
            % clim([0,0.5]);
            % colormap jet
            plot(zeros(size(linspace(0,max(colorCovMat(:)),2))), linspace(0,max(colorCovMat(:)),2),'k')
            xlabel(sprintf("Fold Change in FR (%s-%s)/%s",stateNames{s2},stateNames{s1},stateNames{s1}))
            ylabel(sprintf("Spatial Coverage"))
            % set(gca, 'XScale', 'log');
        end
    end
end
