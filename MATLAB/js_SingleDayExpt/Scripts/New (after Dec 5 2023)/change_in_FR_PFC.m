% Script to analyze changes in firing rate across epochs. Specifcially,
% comparing PFC participation in ripples to changes in REM firing rates.
% Michael Satchell 01/23/24

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


%% Comparing change in FR during select sleep states.
%  Adapted from code in spatial_tuning_PFC.m 


FRChng = [];
states = ["SWS","REM","ripple"]; % must be "SWS","REM","ripple","run","still"
stateEpochs = [5,13]; % Which epochs to calculate change in FR between. Note
% that the difference will be calculated as stateEpoch(2) - stateEpoch(1).

if numel(stateEpochs) > 2
    error("Only 2 epochs are permitted in stateEpochs, as the difference in" + ...
        " FR is being calculated.")
end

% This loop fills FRChng, which holds the change in FR between the two
% stateEpoch epochs for multiple states.
for r = 1:size(C_nrninfo,2)

    for n = 1:size(C_nrninfo{1,r},1)

        nrn = C_nrninfo{1,r}{n,1};
        
        % This conditional can be changed to decide which PFC neurons get
        % analyzed.
        if strcmp(nrn.area,"PFC") && strcmp(nrn.type,"Pyr") && strcmp(nrn.LMRVtype,"beh-LMRV")

            if all([nrn.eHasSpikeData(1,stateEpochs(1)), ...
                    nrn.eHasSpikeData(1,stateEpochs(2))])
            
                stateEpochFRs = NaN(numel(stateEpochs),numel(states));
                % Loop through epoch indices
                for ei = 1:numel(stateEpochs) % Should only loop twice.

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


% s1 and s2 are indices into 'states' variable
s1 = 3;
s2 = 2;
nanIdxs = logical(isnan(FRChng(:,s1))+isnan(FRChng(:,s2)));
s1FRs = FRChng(~nanIdxs,s1); s2FRs = FRChng(~nanIdxs,s2);
[ccoeffs, pval] = corrcoef(s1FRs,s2FRs);

figure;
hold on
xline(0,'--k')
yline(0,"--k")
plot(s1FRs,s2FRs,'o',Color=[144, 20, 222]/222) % purple
title(sprintf("beh-LMRV Change in FR: %s vs %s",states(s2),states(s1)))
ylabel(sprintf("Change in %s FR (Epoch %d - %d)",states(s2),stateEpochs(2),stateEpochs(1)))
xlabel(sprintf("Change in %s FR (Epoch %d - %d)",states(s1),stateEpochs(2),stateEpochs(1)))

fitL = lsline;
lstxt = sprintf("corr=%.1d \n p=%.1d",ccoeffs(2),pval(2));
text(0,fitL.YData(2),lstxt)



%% Compare average FRs for multiple epochs. Plots similar graphs to the 
% ones in the preceding section, but stores FRs differently so that they
% can be averaged across multiple epochs.
% Adapted from spatial_tuning_PFC.m

states = ["SWS","REM","ripple"]; % must be "SWS","REM","ripple","run","still"
epochs = [3,5,7,13,15,17]; % Epochs to store state FR data from. Epochs here must have the states in 'states'.
stateEpochFRs = cell(1,numel(states)); % To hold the FR for each neuron in each epoch

for r = 1:size(C_nrninfo,2)

    for n = 1:size(C_nrninfo{1,r},1)

        nrn = C_nrninfo{1,r}{n,1};
        
        if strcmp(nrn.area,"PFC") && strcmp(nrn.type,"Pyr") && strcmp(nrn.LMRVtype,"beh-LMRV")

            % Loop through states
            for s = 1:numel(states)
                epochFRs = NaN(1,numel(epochs));
                % Find the number of spikes within each state
                for ei = 1:numel(epochs)
                    stateNames = C_combstates{1,r}{1,ei}.stateNames; % state names that exist 
                    % within this epoch
                    stateidx = find(strcmp(stateNames,states(s)));
                    occs = C_combstates{1,r}{1,epochs(ei)}.sepData{stateidx,1};
                    
                    if ~isempty(occs) && nrn.eHasSpikeData(epochs(ei))
                        spikes = nrn.eSpikeData{1,epochs(ei)}(:,1);
    
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


% s1 and s2 are indices into 'states' variable
s1 = 3; s2 = 2;
% Epoch indices for grouping and averaging. Note that change in FR is
% calculated as mean(group2) - mean(group1).
eGroup1 = [1,2,3]; eGroup2 = [4,5,6];

s1FRs = stateEpochFRs{1,s1}; 
s2FRs = stateEpochFRs{1,s2};

% Average across epochs that are grouped together for each state
s1MeanChng = mean(s1FRs(:,eGroup2),2,'omitnan') - mean(s1FRs(:,eGroup1),2,'omitnan');
s2MeanChng = mean(s2FRs(:,eGroup2),2,'omitnan') - mean(s2FRs(:,eGroup1),2,'omitnan');

% Remove any NaNs
nanIdxs = logical(isnan(s1MeanChng)+isnan(s2MeanChng));
s1MeanChng = s1MeanChng(~nanIdxs); 
s2MeanChng = s2MeanChng(~nanIdxs); 

[ccoeffs, pval] = corrcoef(s1MeanChng,s2MeanChng);

figure;
hold on
xline(0,'--k')
yline(0,"--k")
plot(s1MeanChng,s2MeanChng,'o',Color=[222, 144, 20]/222) % orange
%title(sprintf("beh-LMRV Change in FR: %s vs %s",states(s2),states(s1)))
title(sprintf("beh-LMRV Change in FR: %s vs %s",states(s2),states(s1)))
ylabel(sprintf("Change in %s FR ( mean(%s) - mean(%s) )",states(s2),...
    join(string(epochs(eGroup2))), join(string(epochs(eGroup1)))))
xlabel(sprintf("Change in %s FR ( mean(%s) - mean(%s) )",states(s1),...
    join(string(epochs(eGroup2))), join(string(epochs(eGroup1)))))
fitL = lsline;
lstxt = sprintf("corr=%.1d \n p=%.1d",ccoeffs(2),pval(2));
text(0,fitL.YData(2),lstxt)
