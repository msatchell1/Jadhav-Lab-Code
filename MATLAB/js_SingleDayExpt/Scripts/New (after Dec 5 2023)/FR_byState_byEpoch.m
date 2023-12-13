% Script to analyze firing rate data of neurons in different states.
% Additionally this will be split by epoch to look for progressions over
% learning.
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




%% Split spikes up depending on the state they occur during


% Assmues all rats have same number of epochs
for e = 1:size(C_combstates{1,1},2)
    
    % All rats and epochs have the same stateNames
    stateNames = C_combstates{1,1}{1,e}.stateNames;

    % % The index number into stateNames of which states exist for this
    % % epoch. Assuming identical across rats.
    % stateIdxs = find(~cellfun(@isempty,C_combstates{1,1}{1,e}.sepData));

    % Because some rats will have different states (some rats have no REM
    % in epoch 1), I can't just use rat 1 as a reference for what states
    % exist. Instead I just include all possible states in each epoch:
    if mod(e,2) % rest epochs
        stateIdxs = [1:numel(stateNames)]';
    else % beh epochs
        stateIdxs = [3;4;5]; % ripple, run, and still respectively. 
    end

    % Create an array to keep track of which state combinations have
    % been plotted. The number unique pairs is n*(n-1)/2
    numUnqPrs = numel(stateIdxs)*(numel(stateIdxs)-1)/2;
    stateHist = [NaN,NaN]; % Can't be empty to use ismember

    tl = tiledlayout(ceil(numUnqPrs/ceil(sqrt(numUnqPrs))),ceil(sqrt(numUnqPrs)));
    title(tl, sprintf("Epoch %i all Rats",e))
    tlaxes = [];

    for s1 = stateIdxs'
        % Loop through a second set of states to get every combination
        for s2 = stateIdxs'
            % Exclude comparing state with self or duplicating
            % comparisons.
            if s1 ~= s2 && ~any(ismember(stateHist,[s1,s2 ; s2,s1],'rows'))
                stateHist = [stateHist; s1,s2]; % record state comparison made.

                tlaxes(end+1) = nexttile;
                hold on
                xlabel(sprintf("FR in %s (Hz)",stateNames{s1}))
                ylabel(sprintf("FR in %s (Hz)",stateNames{s2}))
                set(gca, 'XScale', 'log', 'YScale', 'log');
                plot([0.01,40],[0.01,40],'k',HandleVisibility='off')
                pltcolors = {[152,8,230]/230, [8,52,230]/230};

                % It's important to note that neurons in nrninfo_MES
                % without an area field also have no spike data in any
                % epochs. There are actually quite a few cells like this in
                % every rat.
                FRs = [];
                areas = [];
                for r = 1:size(C_nrninfo,2)
                    % Get the state occurance times
                    s1Occs = C_combstates{1,r}{1,e}.sepData{s1,1};
                    s2Occs = C_combstates{1,r}{1,e}.sepData{s2,1};
                    if ~any([isempty(s1Occs),isempty(s2Occs)]) % Makes sure both
                        % states have at least one occurance for this
                        % rat/epoch.
                        s1s2Occs = [s1Occs; s2Occs]; % Combine so I can loop easier later.
                        % Loop through neurons and calculate FRs.
                        for n = 1:size(C_nrninfo{1,r},1)
                        
                            S_nrn = C_nrninfo{1,r}{n,1};
                            if S_nrn.eHasSpikeData(1,e) % If nrn has spike data for this epoch 
                                % if strcmp(S_nrn.type,"Int")
                                if strcmp(S_nrn.area,"PFC")
                                    % nrnColor = pltcolors{1}; %"red";
                                    areas = [areas; "PFC"];
                                % elseif strcmp(S_nrn.type,"Pyr")
                                elseif strcmp(S_nrn.area,"CA1")
                                    % nrnColor = pltcolors{2}; %"blue";
                                    areas = [areas; "CA1"];
                                end
        
                                spikes = S_nrn.eSpikeData{1,e}(:,1);
                                % Holds the number of spikes for a given neuron
                                % that occur during each state occurrance.
                                spikesInOcc = zeros(size(s1s2Occs,1),1);
        
                                % Loop through occurances of the two states s1 and
                                % s2, counting the number of spikes that fall into
                                % each occurrance during this epoch.
                                for o = 1:size(s1s2Occs,1)
                                    isInOcc = (s1s2Occs(o,1) <= spikes) & (s1s2Occs(o,2) > spikes);
                                    spikesInOcc(o,1) = sum(isInOcc);
                                end
        
                                % Sum respective spikes for each state
                                s1Sum = sum(spikesInOcc(s1s2Occs(:,3)==s1));
                                s2Sum = sum(spikesInOcc(s1s2Occs(:,3)==s2));
                                % total durations of states
                                s1Dur = sum(s1Occs(:,2) - s1Occs(:,1));
                                s2Dur = sum(s2Occs(:,2) - s2Occs(:,1));
        
                                % store FRs in matrix
                                FRs = [FRs; s1Sum/s1Dur,s2Sum/s2Dur];
                            end
                        end
                    end
                end
                scatter(FRs(areas=="CA1",1),FRs(areas=="CA1",2),[],pltcolors{2},'.')
                scatter(FRs(areas=="PFC",1),FRs(areas=="PFC",2),[],pltcolors{1},'.')
                if size(stateHist,1) == 2 % First plot only
                legend({"CA1","PFC"},Location='best')
                end
            end
        end
    end
    
    linkaxes(tlaxes,'xy')
    pause
end


%% Look at only PFC cells and color code by LMRV


for e = 1:size(C_combstates{1,1},2)
    
    % All rats and epochs have the same stateNames
    stateNames = C_combstates{1,1}{1,e}.stateNames;

    % % The index number into stateNames of which states exist for this
    % % epoch. Assuming identical across rats.
    % stateIdxs = find(~cellfun(@isempty,C_combstates{1,1}{1,e}.sepData));

    % Because some rats will have different states (some rats have no REM
    % in epoch 1), I can't just use rat 1 as a reference for what states
    % exist. Instead I just include all possible states in each epoch:
    if mod(e,2) % rest epochs
        stateIdxs = [1:numel(stateNames)]';
    else % beh epochs
        stateIdxs = [3;4;5]; % ripple, run, and still respectively. 
    end

    % Create an array to keep track of which state combinations have
    % been plotted. The number unique pairs is n*(n-1)/2
    numUnqPrs = numel(stateIdxs)*(numel(stateIdxs)-1)/2;
    stateHist = [NaN,NaN]; % Can't be empty to use ismember
    
    figure;
    tl = tiledlayout(ceil(numUnqPrs/ceil(sqrt(numUnqPrs))),ceil(sqrt(numUnqPrs)));
    title(tl, sprintf("Epoch %i all Rats PFC",e))
    tlaxes = [];

    for s1 = stateIdxs'
        % Loop through a second set of states to get every combination
        for s2 = stateIdxs'
            % Exclude comparing state with self or duplicating
            % comparisons.
            if s1 ~= s2 && ~any(ismember(stateHist,[s1,s2 ; s2,s1],'rows'))
                stateHist = [stateHist; s1,s2]; % record state comparison made.

                tlaxes(end+1) = nexttile;
                hold on
                xlabel(sprintf("FR in %s (Hz)",stateNames{s1}))
                ylabel(sprintf("FR in %s (Hz)",stateNames{s2}))
                % set(gca, 'XScale', 'log', 'YScale', 'log');
                pltcolors = {[152,8,230]/230, [230,10,2]/230, [27,230,230]/230};

                % It's important to note that neurons in nrninfo_MES
                % without an area field also have no spike data in any
                % epochs. There are actually quite a few cells like this in
                % every rat.
                FRs = [];
                LMRV = [];
                isInh = [];
                pkEpoch = [];
                for r = 1:size(C_nrninfo,2)
                    % Get the state occurance times
                    s1Occs = C_combstates{1,r}{1,e}.sepData{s1,1};
                    s2Occs = C_combstates{1,r}{1,e}.sepData{s2,1};
                    if ~any([isempty(s1Occs),isempty(s2Occs)]) % Makes sure both
                        % states have at least one occurance for this
                        % rat/epoch.
                        s1s2Occs = [s1Occs; s2Occs]; % Combine so I can loop easier later.
                        % Loop through neurons and calculate FRs.
                        for n = 1:size(C_nrninfo{1,r},1)
                        
                            S_nrn = C_nrninfo{1,r}{n,1};
                            if S_nrn.eHasSpikeData(1,e) && strcmp(S_nrn.area,"PFC") % If nrn has spike data for this epoch 
                                % if strcmp(S_nrn.type,"Int")
                                if strcmp(S_nrn.LMRVtype,"rest-LMRV")
                                    LMRV = [LMRV; 1];
                                % elseif strcmp(S_nrn.type,"Pyr")
                                elseif strcmp(S_nrn.LMRVtype,"beh-LMRV")
                                    LMRV = [LMRV; 2];
                                else
                                    LMRV = [LMRV; 0];
                                end

                                % Find interneurons
                                if strcmp(S_nrn.type,"Int")
                                    isInh = [isInh; 1];
                                else
                                    isInh = [isInh; 0];
                                end

                                % Find epoch with highest FR
                                [mx,mxI] = max(S_nrn.eFR);
                                pkEpoch = [pkEpoch; mxI];
        
                                spikes = S_nrn.eSpikeData{1,e}(:,1);
                                % Holds the number of spikes for a given neuron
                                % that occur during each state occurrance.
                                spikesInOcc = zeros(size(s1s2Occs,1),1);
        
                                % Loop through occurances of the two states s1 and
                                % s2, counting the number of spikes that fall into
                                % each occurrance during this epoch.
                                for o = 1:size(s1s2Occs,1)
                                    isInOcc = (s1s2Occs(o,1) <= spikes) & (s1s2Occs(o,2) > spikes);
                                    spikesInOcc(o,1) = sum(isInOcc);
                                end
        
                                % Sum respective spikes for each state
                                s1Sum = sum(spikesInOcc(s1s2Occs(:,3)==s1));
                                s2Sum = sum(spikesInOcc(s1s2Occs(:,3)==s2));
                                % total durations of states
                                s1Dur = sum(s1Occs(:,2) - s1Occs(:,1));
                                s2Dur = sum(s2Occs(:,2) - s2Occs(:,1));
        
                                % store FRs in matrix
                                FRs = [FRs; s1Sum/s1Dur,s2Sum/s2Dur];
                            end
                        end
                    end
                end
                plot([0.01,max(FRs)],[0.01,max(FRs)],'k',HandleVisibility='off')
                scatter(FRs(LMRV==0,1),FRs(LMRV==0,2),30,pltcolors{1},'.')
                scatter(FRs(LMRV==1,1),FRs(LMRV==1,2),60,pltcolors{2},'.')
                scatter(FRs(LMRV==2,1),FRs(LMRV==2,2),60,pltcolors{3},'.')
                text(FRs(isInh==1,1),FRs(isInh==1,2),"In")
                text(FRs(LMRV==2,1),FRs(LMRV==2,2),num2str(pkEpoch(LMRV==2)),HorizontalAlignment="right")
                if size(stateHist,1) == 2 % First plot only
                legend({"non-LMRV","rest-LMRV","beh-LMRV"},Location='best')
                end
            end
        end
    end
    
    linkaxes(tlaxes,'xy')
    pause
end





%% Determine task-involved cells and look at their FR during NREM and REM
% Could a way to determine task-involved cells be based simply on FR during
% running on track?


%% I could just look at PFC cells that participste in ripples consistently
% and see if they have the biggest increase in REM FR by the final epochs.


%% Measure the difference in FR between REM and SWS for beh-LMRV and all PFC 
% cells for each epoch to try and quantify the differences I see in the FR
% scatter plots. I will compare the distance from the diagonal to the
% firing rate during REM and look for differences in the beh-LMRV cells.

% state indices in C_combstates: {SWS; REM; ripple; run; still}
s1 = 1; 
s2 = 2;

for e = 1:size(C_combstates{1,1},2)
    % Plots only proper epochs for the given states
    if (mod(e,2) && any(ismember([1,2,3],[s1,s2]))) || ...
           (~mod(e,2) && ~any(ismember([1,2,3],[s1,s2])))
        
        % All rats and epochs have the same stateNames
        stateNames = C_combstates{1,1}{1,e}.stateNames;
    
        % % The index number into stateNames of which states exist for this
        % % epoch. Assuming identical across rats.
        % stateIdxs = find(~cellfun(@isempty,C_combstates{1,1}{1,e}.sepData));
        
        figure;
        title(sprintf("Epoch %i PFC Pyr | %s vs %s",e,stateNames{s1},stateNames{s2}))
    
        hold on
        xlabel(sprintf("FR in %s (Hz)",stateNames{s2}))
        ylabel(sprintf("FR Preference towards %s (Hz)",stateNames{s2}))
        % set(gca, 'XScale', 'log', 'YScale', 'log');
        pltcolors = {[152,8,230]/230, [230,10,2]/230, [27,230,230]/230};
    
        % It's important to note that neurons in nrninfo_MES
        % without an area field also have no spike data in any
        % epochs. There are actually quite a few cells like this in
        % every rat.
        FRs = [];
        LMRV = [];
        pkEpoch = [];
        distDiag = []; % Distance from the diagonal of equal firing rate in both states.
    
        for r = 1:size(C_nrninfo,2)
            % Get the state occurance times
            s1Occs = C_combstates{1,r}{1,e}.sepData{s1,1};
            s2Occs = C_combstates{1,r}{1,e}.sepData{s2,1};
            if ~any([isempty(s1Occs),isempty(s2Occs)]) % Makes sure both
                % states have at least one occurance for this
                % rat/epoch.
                s1s2Occs = [s1Occs; s2Occs]; % Combine so I can loop easier later.
                % Loop through neurons and calculate FRs.
                for n = 1:size(C_nrninfo{1,r},1)
                
                    S_nrn = C_nrninfo{1,r}{n,1};
                    if S_nrn.eHasSpikeData(1,e) && strcmp(S_nrn.area,"PFC") ...
                            && ~strcmp(S_nrn.type,"Int") % If nrn has spike data for this epoch 
                        % if strcmp(S_nrn.type,"Int")
                        if strcmp(S_nrn.LMRVtype,"rest-LMRV")
                            LMRV = [LMRV; 1];
                        % elseif strcmp(S_nrn.type,"Pyr")
                        elseif strcmp(S_nrn.LMRVtype,"beh-LMRV")
                            LMRV = [LMRV; 2];
                        else
                            LMRV = [LMRV; 0];
                        end
    
                        % Find epoch with highest FR
                        [mx,mxI] = max(S_nrn.eFR);
                        pkEpoch = [pkEpoch; mxI];
    
                        spikes = S_nrn.eSpikeData{1,e}(:,1);
                        % Holds the number of spikes for a given neuron
                        % that occur during each state occurrance.
                        spikesInOcc = zeros(size(s1s2Occs,1),1);
    
                        % Loop through occurances of the two states s1 and
                        % s2, counting the number of spikes that fall into
                        % each occurrance during this epoch.
                        for o = 1:size(s1s2Occs,1)
                            isInOcc = (s1s2Occs(o,1) <= spikes) & (s1s2Occs(o,2) > spikes);
                            spikesInOcc(o,1) = sum(isInOcc);
                        end
    
                        % Sum respective spikes for each state
                        s1Sum = sum(spikesInOcc(s1s2Occs(:,3)==s1));
                        s2Sum = sum(spikesInOcc(s1s2Occs(:,3)==s2));
                        % total durations of states
                        s1Dur = sum(s1Occs(:,2) - s1Occs(:,1));
                        s2Dur = sum(s2Occs(:,2) - s2Occs(:,1));
                        
                        % store FRs in matrix
                        s1FR = s1Sum/s1Dur;
                        s2FR = s2Sum/s2Dur;
                        FRs = [FRs; s1FR,s2FR];
    
                        % Calculate the distance of each point in FR space from
                        % the diagonal line representing equal FR in both
                        % states. Note this distance equals half the distance
                        % beween the point and its reflection.
                        d = sqrt((s1FR-s2FR)^2 + (s2FR-s1FR)^2)/2;
                        % But d is always positive and we want to define the
                        % s2 direction as positive and s1 as negative,
                        % making the diagonal the new x-axis.
                        if s2FR < s1FR
                            d = -d;
                        end
                        distDiag = [distDiag; d];
    
                        
                    end
                end
            end
        end
        % Sort the firing rates in state 2
        [sorts2FRs,sorts2Idx] = sort(FRs(:,2));
        % Sort the distances to match
        sortDist = distDiag(sorts2Idx);
        plot(sorts2FRs(LMRV==0),sortDist(LMRV==0),Color=pltcolors{1})
        plot(sorts2FRs(LMRV==2),sortDist(LMRV==2),Color=pltcolors{3})
    
        legend({"non-LMRV","beh-LMRV"},Location='best')
    
        pause
    end
end
% I need to make this a histogram counting the number of cells in each FR
% bin.

%% NOTE: remove ripples from SWS and still states (check if they are in still)




