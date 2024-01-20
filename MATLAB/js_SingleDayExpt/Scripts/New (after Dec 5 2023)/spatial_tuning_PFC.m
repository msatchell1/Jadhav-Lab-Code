% Script to asses spatial properties of PFC cells across epochs.
% Michael Satchell 11/02/23


%% Load data

clearvars;

dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
loadRats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01',
% 'spikes01','tetinfo','linfields01','rippletime01','pos01','nrninfo_MES','combstates_MES'
filetypes = {'nrninfo_MES','pos01','behavperform','combstates_MES'};

C_alldata = load_data(dataDir,loadRats,filetypes); % Load data to struct



if ~any(contains(filetypes,'pos01'))
    error("pos01 data must be loaded to run this analysis.")
end
if ~any(contains(filetypes,'behavperform'))
    error("behavperform data must be loaded to run this analysis.")
end
if ~any(contains(filetypes,'nrninfo_MES'))
    error("nrninfo_MES data must be loaded to run this analysis.")
end
if ~any(contains(filetypes,'combstates_MES'))
    error("combstates_MES data must be loaded to run this analysis.")
end

C_behper = C_alldata(find(contains(filetypes,'behavperform')),:);
C_nrninfo = C_alldata(find(contains(filetypes,'nrninfo_MES')),:);
C_pos = C_alldata(find(contains(filetypes,'pos01')),:);
C_combstates = C_alldata(find(contains(filetypes,'combstates_MES')),:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;
brainAreas = {'CA1','PFC'};


%% The portion of an entire trajectory that a cell has a firing field over is
% the spatial coverage of that cell. Spatial coverage differs, naturally,
% between trajectories for the same cell. I want to determine if there is
% much variation between trajectories or if cells have similar spatial
% coverage on all trajectories. This will tell me if I can combine
% trajectory coverages into one number and represent that cell accurately.

% To do this, I will make a histogram of the difference between the largest
% and smallest coverage values for each epoch of each PFC pyr neuron.
spatCovSpan = []; % span of spatial coverage for each epoch (highest-lowest traj val)
spatCovAvg = []; % average value across trajs for each epoch
spatCovChng = []; % difference between first and last epoch in avg coverage value
for r = 1:size(C_nrninfo,2)
    
    for n = 1:size(C_nrninfo{1,r},1)

        nrn = C_nrninfo{1,r}{n,1};
        
        if strcmp(nrn.area,"PFC") && strcmp(nrn.type,"Pyr") && strcmp(nrn.LMRVtype,"non-LMRV")
            covs = nrn.eTrajCoverage;
            isCovExist = cellfun(@(x) ~isempty(x), covs); % Epochs that have spat cov vals
            if sum(isCovExist) > 1 % If at least two epochs have spatial coverage values
                covExist = covs(isCovExist);
                avgCov = cellfun(@(x) mean(x,2,'omitnan'), covExist);
                spatCovAvg = [spatCovAvg, avgCov];
                spatCovChng(end+1) = avgCov(end) - avgCov(1);
                for e = 1:size(covExist,2)
                    spatCovSpan(end+1) = max(covExist{1,e}) - min(covExist{1,e});
                end
            end
        end


    end


end

figure
histogram(spatCovSpan,20)
title("Spatial Coverage Span per Epoch all PFC Pyr non-LMRV")
xlabel("max(cov) - min(cov)")
ylabel("count")

figure
histogram(spatCovAvg,20)
title("Avg Spatial Coverage all Epochs all PFC Pyr non-LMRV")
xlabel("mean cov across trajs")
ylabel("count")

figure
histogram(spatCovChng,30)
title("Change in Spatial Coverage all PFC Pyr non-LMRV")
xlabel("(mean last epoch) - (mean first epoch)")
ylabel("count")


%% Comparing change in sptaial coverage to change in FR during select sleep
% states.
%  think it might make sense to compare with the second, because this is 
% after the initial effects of the track, but before sleep consolidation 
% effects have done their part. Thus what I am comparing is the effect of 
% learning/consolidation, and not novelty.  

% Note this analysis only considers cells that have both more than 1 epoch
% with spatial coverage values and has spike data in the epochs the FR is
% being calculated between.


spatCovChng = []; % difference between first and last epoch in avg coverage value
FRChng = [];
states = ["SWS","REM","ripple"]; % must be "SWS","REM","ripple","run","still"
stateEpochs = [3,9]; % Which epochs to calculate change in FR between. Note
% that the difference will be calculated as stateEpoch(2) - stateEpoch(1).

if numel(stateEpochs) > 2
    error("Only 2 epochs are permitted in stateEpochs, as the difference in" + ...
        " FR is being calculated.")
end


% This loo fills two matrices: spatCovChng which holds the
% change in spatial coverage for each neuron between the first and last
% available epochs, and FRChng which holds the change in FR between the two
% stateEpoch epochs for multiple states.
for r = 1:size(C_nrninfo,2)

    ratStateFRs = cell(numel(states),numel(stateEpochs)); % cell matrix to hold 
    % the occurance times for each state and each epoch.
        
    stateNames = C_combstates{1,r}{1,e}.stateNames; % state names that exist 
    % within this epoch

    % Next step: grab occurances from each epoch and state and store in
    % stateOccs. Then go into the neuron loop and determine number of
    % spikes in each state and epoch. Divide by total epoch time to get
    % FR in that state for each epoch. Then take the difference of this
    % FR between the two epochs to get a single number per neuron.
    % Store this, along with the change in spatial coverage, then plot
    % them against each other and look for correlations. Then repeat
    % with LMRV and beh-LMRV neurons only.

    for n = 1:size(C_nrninfo{1,r},1)

        nrn = C_nrninfo{1,r}{n,1};
        
        if strcmp(nrn.area,"PFC") && strcmp(nrn.type,"Pyr") && ~strcmp(nrn.LMRVtype,"non-LMRV")
            covs = nrn.eTrajCoverage;
            isCovExist = cellfun(@(x) ~isempty(x), covs); % Epochs that have spat cov vals
            
            % If at least two epochs have spatial coverage values and the
            % proper epochs have spike data
            if sum(isCovExist) > 1 && all([nrn.eHasSpikeData(1,stateEpochs(1)), ...
                    nrn.eHasSpikeData(1,stateEpochs(2))])
                
                covExist = covs(isCovExist);
                avgCov = cellfun(@(x) mean(x,2,'omitnan'), covExist);
                spatCovChng = [spatCovChng; avgCov(end) - avgCov(1)];
            
                stateEpochFRs = NaN(numel(stateEpochs),numel(states));
                % Loop through epoch indices
                for ei = 1:numel(stateEpochs)

                    % Find the number of spikes within each state
                    for s = 1:numel(states)
                        stateidx = find(strcmp(stateNames,states(s)));
                        occs = C_combstates{1,r}{1,stateEpochs(ei)}.sepData{stateidx,1};
                        
                        if ~isempty(occs)
                            spikes = nrn.eSpikeData{1,stateEpochs(ei)}(:,1);
        
                            % Holds the number of spikes for a given neuron
                            % that occur during each state occurrance.
                            numSpikesInOcc = zeros(size(occs,1),1);
            
                            % Loop through occurances of the state,
                            % counting the number of spikes that fall into
                            % each occurrance.
                            for o = 1:size(occs,1)
                                % spike times greater than state start time and less
                                % than state end time.
                                isInOcc = (occs(o,1) <= spikes) & (occs(o,2) > spikes);
                                numSpikesInOcc(o,1) = sum(isInOcc);
                            end
            
                            spikeSum = sum(numSpikesInOcc);
                            stateDur = sum(occs(:,2) - occs(:,1)); % total duration of state
                            stateEpochFRs(ei,s) = spikeSum/stateDur; % store FR
                        end
                    end
                end

                FRChng = [FRChng; stateEpochFRs(2,:)-stateEpochFRs(1,:)];
            end
        end
    end
end


for s = 1:numel(states)
    figure;
    plot(spatCovChng,FRChng(:,s),'o')
    title(sprintf("Change in %s FR vs Change in Spatial Coverage",states(s)))
    ylabel(sprintf("Change in FR (Epoch %d - %d)",stateEpochs(2),stateEpochs(1)))
    xlabel("Change in spatial coverage")
end