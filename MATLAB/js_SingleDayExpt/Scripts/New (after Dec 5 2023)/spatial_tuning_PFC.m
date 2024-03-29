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

% As a reminder: There are four columns here representing the four W-track trajectories.
% Col 1: out right, col 2: in right, col 3: out left, col 4: in left.

% To do this, I will make a histogram of the difference between the largest
% and smallest coverage values for each epoch of each PFC pyr neuron.
spatCovSpan = []; % span of spatial coverage for each epoch (highest-lowest traj val)
spatCovAvg = []; % average value across trajs for each epoch
spatCovChng = []; % difference between first and last epoch in avg coverage value
cellType = "non-LMRV";
trajFR = []; % Firing rates on all trajectories for all epochs and cells
trajSC = [];
for r = 1:size(C_nrninfo,2)
    
    for n = 1:size(C_nrninfo{1,r},1)

        nrn = C_nrninfo{1,r}{n,1};
        
        if strcmp(nrn.area,"PFC") && strcmp(nrn.type,"Pyr") && strcmp(nrn.LMRVtype,cellType)
            covs = nrn.eTrajCoverage;
            isCovExist = cellfun(@(x) ~isempty(x), covs); % Epochs that have spat cov vals
            if sum(isCovExist) > 1 % If at least two epochs have spatial coverage values
                covExist = covs(isCovExist);
                avgCov = cellfun(@(x) mean(x,2,'omitnan'), covExist);
                spatCovAvg = [spatCovAvg, avgCov];
                spatCovChng(end+1) = avgCov(end) - avgCov(1);
                for e = 1:size(covExist,2)
                    spatCovSpan(end+1) = max(covExist{1,e}) - min(covExist{1,e});
                    trajFR = [trajFR; nrn.eTrajFR{1,e}];
                end
            end
        end


    end


end

% figure
% histogram(spatCovSpan,20)
% title(sprintf("Spatial Coverage Span per Epoch all PFC Pyr %s",LMRVtype))
% xlabel("max(cov) - min(cov)")
% ylabel("count")
% 
% figure
% histogram(spatCovAvg,20)
% title(sprintf("Avg Spatial Coverage all Epochs all PFC Pyr %s",LMRVtype))
% xlabel("mean cov across trajs")
% ylabel("count")
% 
% figure
% histogram(spatCovChng,30)
% title(sprintf("Change in Spatial Coverage all PFC Pyr %s",LMRVtype))
% xlabel("(mean last epoch) - (mean first epoch)")
% ylabel("count")

figure
histogram(trajFR(:),50)
title(sprintf("Trajectory FR all PFC Pyr all Epochs %s",cellType))
xlabel("Rate (Hz)")
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
stateEpochs = [3,17]; % Which epochs to calculate change in FR between. Note
% that the difference will be calculated as stateEpoch(2) - stateEpoch(1).

if numel(stateEpochs) > 2
    error("Only 2 epochs are permitted in stateEpochs, as the difference in" + ...
        " FR is being calculated.")
end

%allTrajCov = [];

% This loo fills two matrices: spatCovChng which holds the
% change in spatial coverage for each neuron between the first and last
% available epochs, and FRChng which holds the change in FR between the two
% stateEpoch epochs for multiple states.
for r = 1:size(C_nrninfo,2)

    for n = 1:size(C_nrninfo{1,r},1)

        nrn = C_nrninfo{1,r}{n,1};
        
        % This conditional can be changed to decide which PFC neurons get
        % analyzed.
        if strcmp(nrn.area,"PFC") && strcmp(nrn.type,"Pyr") && strcmp(nrn.LMRVtype,"beh-LMRV")
            covs = nrn.eTrajCoverage;
            isCovExist = cellfun(@(x) ~isempty(x), covs); % Epochs that have spat cov vals
            %allTrajCov = [allTrajCov; nrn.eTrajCoverage{1,2},NaN,nrn.eTrajCoverage{1,16}];
            % If at least two epochs have spatial coverage values and the
            % proper epochs have spike data
            if sum(isCovExist) > 1 && all([nrn.eHasSpikeData(1,stateEpochs(1)), ...
                    nrn.eHasSpikeData(1,stateEpochs(2))])
                %[find(isCovExist,1,"first"),find(isCovExist,1,"last")]
                covExist = covs(isCovExist);
                %avgCov = cellfun(@(x) mean(x,2,'omitnan'), covExist);

                % I want to find the trajectory that has the greatest drop
                % in coverage from the first and last beh epoch
                maxCovChng = min(covExist{end}-covExist{1});
                %spatCovChng = [spatCovChng; avgCov(end) - avgCov(1)];
                spatCovChng = [spatCovChng; maxCovChng];
            
                stateEpochFRs = NaN(numel(stateEpochs),numel(states));
                % Loop through epoch indices
                for ei = 1:numel(stateEpochs)

                    stateNames = C_combstates{1,r}{1,ei}.stateNames; % state names that exist 
                    % within this epoch 

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
    hold on
    xline(0,'--k')
    yline(0,"--k")
    plot(spatCovChng,FRChng(:,s),'o')
    title(sprintf("Change in %s FR vs Change in Spatial Coverage",states(s)))
    ylabel(sprintf("Change in FR (Epoch %d - %d)",stateEpochs(2),stateEpochs(1)))
    xlabel("Largest drop in coverage")
end


%% Compare average FR for multiple epochs to spatial coverage change


spatCovChng = []; % difference between first and last epoch in avg coverage value.
states = ["SWS","REM","ripple"]; % must be "SWS","REM","ripple","run","still"
epochsToAvg = [17]; % Different from stateEpochs above. Instead, all these epochs
% will have the firing rates of each state calculated within them, and then
% be averaged together. Epochs here must have the states in 'states'.
stateEpochFRs = cell(1,numel(states)); % To hold the FR for each neuron in each epoch

for r = 1:size(C_nrninfo,2)

    for n = 1:size(C_nrninfo{1,r},1)

        nrn = C_nrninfo{1,r}{n,1};
        
        if strcmp(nrn.area,"PFC") && strcmp(nrn.type,"Pyr") && strcmp(nrn.LMRVtype,"beh-LMRV")
            covs = nrn.eTrajCoverage;
            isCovExist = cellfun(@(x) ~isempty(x), covs); % Epochs that have spat cov vals
            % If at least two epochs have spatial coverage values and the
            % proper epochs have spike data
            if sum(isCovExist) > 1

                covExist = covs(isCovExist);
                avgCov = cellfun(@(x) mean(x,2,'omitnan'), covExist);

                % I want to find the trajectory that has the greatest drop
                % in coverage from the first and last beh epoch
                maxCovChng = min(covExist{end}-covExist{1});
                spatCovChng = [spatCovChng; avgCov(end) - avgCov(1)];
                % spatCovChng = [spatCovChng; maxCovChng];

                % Loop through epoch indices (not epoch numbers)
                for s = 1:numel(states)
                    epochFRs = NaN(1,numel(epochsToAvg));
                    % Find the number of spikes within each state
                    for ei = 1:numel(epochsToAvg)
                        stateNames = C_combstates{1,r}{1,ei}.stateNames; % state names that exist 
                        % within this epoch
                        stateidx = find(strcmp(stateNames,states(s)));
                        occs = C_combstates{1,r}{1,epochsToAvg(ei)}.sepData{stateidx,1};
                        
                        if ~isempty(occs) && nrn.eHasSpikeData(epochsToAvg(ei))
                            spikes = nrn.eSpikeData{1,epochsToAvg(ei)}(:,1);
        
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
                            epochFRs(1,ei) = spikeSum/stateDur; % store FR
                        end
                    end
                    stateEpochFRs{1,s} = [stateEpochFRs{1,s}; epochFRs];
                end
            end
        end
    end
end


for s = 1:numel(states)
    avgFR = mean(stateEpochFRs{1,s},2,'omitnan');
    % Sort out neurons that had NaN FRs for all epochs
    nanIdxs = isnan(avgFR);
    avgFR = avgFR(~nanIdxs);
    reducedCovChng = spatCovChng(~nanIdxs); % Important to assign this to a
    % new variable because which neurons have NaN values changes with each
    % state.
    [ccoeffs, pval] = corrcoef(reducedCovChng,avgFR);
    figure;
    hold on
    xline(0,'--k')
    plot(reducedCovChng,avgFR,'o')
    fitL = lsline;
    lstxt = sprintf("corr=%.2d \n p=%.2d",ccoeffs(2),pval(2));
    text((fitL.XData(1)-fitL.XData(2))/2,fitL.YData(2),lstxt)
    title(sprintf("Avg %s FR vs Change in Spatial Coverage",states(s)))
    ylabel("Avg FR over Epochs: "+join(string(epochsToAvg),","))
    xlabel("Avg drop in coverage (all trajs)")
    % xlabel("Largest drop in coverage (single traj)")
end



%% Combining trajectory measures for all neurons separated by epoch
% As a reminder: There are four columns here representing the four W-track trajectories.
% Col 1: out right, col 2: in right, col 3: out left, col 4: in left.

cellType = "CA1";
trajFR = cell(size(C_combstates{1,1})); % Firing rates on all trajectories for all epochs and cells
trajSC = cell(size(C_combstates{1,1})); % Spatial coverage
trajPk = cell(size(C_combstates{1,1})); % Peak firing rate
trajisPC = cell(size(C_combstates{1,1})); % Place cell (1) or not (0)
for r = 1:size(C_nrninfo,2)
    
    for n = 1:size(C_nrninfo{1,r},1)

        nrn = C_nrninfo{1,r}{n,1};
        
        if strcmp(nrn.area,"CA1") && strcmp(nrn.type,"Pyr") %&& strcmp(nrn.LMRVtype,cellType)
    
            for e = 1:size(nrn.eTrajFR,2) % Loops through all epochs
                if ~isempty(nrn.eTrajFR{1,e}) % assume that if this nrn has FR data 
                    % for this epoch, then all other traj measures will
                    % exist.
                    trajFR{1,e} = [trajFR{1,e}; nrn.eTrajFR{1,e}];
                    trajSC{1,e} = [trajSC{1,e}; nrn.eTrajCoverage{1,e}];
                    trajPk{1,e} = [trajPk{1,e}; nrn.eTrajFRPeak{1,e}];
                    trajisPC{1,e} = [trajisPC{1,e}; nrn.eTrajisPC{1,e}];
                end
            end

        end


    end


end

% figure
% tl = tiledlayout(2,1);
% e = 2; % Even epochs only
% numBins = 20;
% 
% nexttile
% histogram(trajSC{1,e}(:,[2,4]),numBins) % inbound trajs
% title("Inbound")
% 
% nexttile
% histogram(trajSC{1,e}(:,[1,3]),numBins) % outbound trajs
% title("Outbound")
% 
% title(tl,sprintf("Spatial Coverage on Epoch %d for %s",e,cellType))
% xlabel(tl,"Spatial Coverage")
% ylabel(tl,"Count")
% linkaxes


% figure
% tl = tiledlayout(2,1);
% e = 2; % Even epochs only
% numBins = 30;
% 
% nexttile
% histogram(trajFR{1,e}(:,[2,4]),numBins) % inbound trajs
% title("Inbound")
% 
% nexttile
% histogram(trajFR{1,e}(:,[1,3]),numBins) % outbound trajs
% title("Outbound")
% 
% title(tl,sprintf("Mean Trajectory FR on Epoch %d for %s",e,cellType))
% xlabel(tl,"Rate (Hz)")
% ylabel(tl,"Count")
% linkaxes



% figure
% tl = tiledlayout(2,1);
% e = 16; % Even epochs only
% dataArray = trajPk{1,e};
% numBins = 30;
% binEdges = linspace(min(dataArray(:)),max(dataArray(:)),numBins);
% 
% nexttile
% histogram(dataArray(:,[2,4]),binEdges) % inbound trajs
% title("Inbound")
% 
% nexttile
% histogram(dataArray(:,[1,3]),binEdges) % outbound trajs
% title("Outbound")
% 
% title(tl,sprintf("Trajectory Peak FR on Epoch %d for %s",e,cellType))
% xlabel(tl,"Rate (Hz)")
% ylabel(tl,"Count")
% linkaxes


 
% figure
% tl = tiledlayout(2,1);
% e = 2; % Even epochs only
% dataArray = trajisPC{1,e};
% 
% binEdges = [-0.5,0.5,1.5,2.5];
% 
% nexttile
% histogram(sum(dataArray(:,[2,4]),2),binEdges) % inbound trajs
% title("Inbound")
% 
% nexttile
% histogram(sum(dataArray(:,[1,3]),2),binEdges) % outbound trajs
% title("Outbound")
% 
% title(tl,sprintf("Number of Place Fields on Epoch %d for %s",e,cellType))
% xlabel(tl,"Number of Place Fields for a Cell")
% ylabel(tl,"Count")
% linkaxes


% % Correlate trajectory peak FR with trajectory spatial coverage
% figure
% e1 = 2;
% trajsToPlot = [2,4];
% data1 = reshape(trajSC{1,e1}(:,trajsToPlot),[],1);
% data2 = reshape(trajPk{1,e1}(:,trajsToPlot),[],1);
% hold on
% % xline(0,'--k')
% % yline(0,"--k")
% plot(data1,data2,'o',Color=[50, 190, 210]/222) % cyan
% title(sprintf("Peak FR and Spatial Coverage for Traj %s, %s",join(string(trajsToPlot)),cellType))
% ylabel(sprintf("Peak Rate (Hz)"))
% xlabel(sprintf("Spatial Coverage"))
% 
% % Plot another epoch on same graph
% e2 = 16;
% data1 = reshape(trajSC{1,e2}(:,trajsToPlot),[],1);
% data2 = reshape(trajPk{1,e2}(:,trajsToPlot),[],1);
% hold on
% plot(data1,data2,'o',Color=[0, 0, 0]/222) % black
% fitL = lsline();
% [fitL.LineWidth] = deal(2,2);
% legend([sprintf("epoch %d",e1),sprintf("epoch %d",e2)])


% figure
% combSC = [];
% for e = 1:size(trajSC,2)
%     combSC = [combSC; trajSC{1,e}];
% end
% medSC = median(combSC,2);
% h = histogram(medSC);
% h.FaceColor = 'r';
% title(sprintf("Spatial Coverage For Each Epoch %s",cellType))
% xlabel("Median Coverage Across All 4 Trajectories")
% ylabel("Count")
% ylim([0,700])
