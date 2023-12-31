function [C_stateFRs,FR_allStates,Occ_allStates,Dur_allStates] = calc_meanrates(brainAreas,C_allStates,C_allSpikes)
%CALC_MEANRATES Calculates the mean firing rates for neurons in the states
%included in C_allstates based on occurances of each state.
% Michael Satchell 09/22/23
%
% NOTE: pyramidal and inhibitory outputs were removed because far too many
% pyramidal cells were being missed, for some reason. I plotted epoch 2
% across all rats and using th isPyrMat on FR_allStates dropped the number
% of cells from 532 to 340, which is far more than should be the case, and
% looking at the epoch 2 FRs, many of the cells dropped were below 7 Hz in
% epoch 2. 
%
% Inputs:
%   brainAreas - cell array of string, brain regions of interest.
%
%   C_allStates - cell array containing loaded state data from each rat.
%   This is used to determine the occurance times for each state.
%
%   C_allSpikes - cell array containing loaded spike data from each rat.
%   This is used to get the mean spike rate data from each neuron.
%
% Outputs:
%   C_stateFRs - Cell array with the same data as FR_allstates, but data
%   from each rat is kept separate.
%
%   C_isPyr - Cell array 1 x (num rats) of logical values
%   indicating if neuron is a pyramidal cell (< 7 Hz mean rate)
%
%   C_isInh - Same as C_isPyr but for inhibitory neurons.
%
%   FR_allStates - (num states) x (num brain areas) cell array. Each cell
%   contains a matrix of concatenated mean rates of size (tot num neurons)
%   x (num epochs). The rows are neurons from C_allspikes concatenated
%   across all rats, which is constant across epochs allowing for a square
%   matrix. NaNs pad the matrix where cells exist either in a different
%   brain area or different epoch.
%
%   isPyrMat - (tot num neurons) x (num epochs) logical matrix indicating
%   if neuron is a pyramidal cell.
%
%   isInhMat - (tot num neurons) x (num epochs) logical matrix indicating
%   if neuron is an inhibitory cell.
%
%   Occ_allStates - (num states) x (num rats) cell array. Each cell
%   contains a cell array of 1 x (num epochs), and within those cells are
%   are the start (col 1) and end (col 2) times for separate occurances of
%   that state.
%
%   Dur_allStates - Matrix of total time duration of states by rat and
%   epoch, i.e. the sum of all occurances.
%   (num states) x  (num rats) x (num epochs).



% Gets epochs that have data for all states. For most states this should be all odd
% epochs for the data I am sorting here, because all the states
% happen during the rest epochs except for rippletimes. I will detect non-empty
% epochs for generality. Note this assumes all rats have the same number of
% epochs as rat 1, which should be true when clip_17_epochs.m has been
% used on C_allstates.
nonemptyEpochs = zeros(size(C_allStates,1),length(C_allStates{1,1}));
for s2 = 1:size(C_allStates,1)
    nonemptyEpochs(s2,:) = ~cellfun(@isempty, C_allStates{s2,1});
end

% To hold firing rates for all states held in C_allstates
FR_allStates = cell(size(C_allStates,1),length(brainAreas));
% Logical matrices indicating pyramidal or inhibitory cell
isPyrMat = [];
isInhMat = [];
flag = 0; % To only add data to the isPyr and isInh mats only once.

% Time duration of states by epoch. (num states) x  (num rats) x (num epochs).
Dur_allStates = zeros(size(C_allStates,1),size(C_allStates,2),size(nonemptyEpochs,2));
% Dur_ratRef = zeros(size(C_allStates,1),size(C_allStates,2),size(nonemptyEpochs,2)); % Array to reference rat number
% All the start and end times of all occurances. (num states) x (num rats).
% Then within each cell are cells for each epoch, and within those is a
% (num occurances) x 2 matrix. Col 1 is starttimes, col 2 endtimes.
Occ_allStates = cell(size(C_allStates,1),size(C_allStates,2));

% Main loop to calculate mean firing rate for all states in all brain areas.
for a = 1:length(brainAreas)
    fprintf("Working on %s... \n",brainAreas{a})
    
    for s2 = 1:size(C_allStates,1) % Loops through all states (e.g. 'sws','rem','ripple'...)
    
        % Gets epochs that have data for this state. This should be identical across rats.
        nonempty = zeros(1,length(C_allStates{1,1}));
    
        for r = 1:size(C_allStates,2) 
            % Indicates which epochs contain data
            nonempty(1,:) = ~cellfun(@isempty, C_allStates{s2,r});
            % To hold mean firing rate data across all epochs for a single rat.
            ratMeans = NaN(length([C_allSpikes{1,r}{1,1}{:}]), length(C_allSpikes{1,r}));
            % To hold 1s and 0s indicating is pyramidal/inhibitory cell
            ratisPyr = zeros(length([C_allSpikes{1,r}{1,1}{:}]), length(C_allSpikes{1,r}));
            ratisInh = zeros(length([C_allSpikes{1,r}{1,1}{:}]), length(C_allSpikes{1,r}));
            % To hold state occurances for all epochs.
            occEpochs = cell(1,length(C_allStates{s2,r}));
    
            for e = 1:length(C_allStates{s2,r})
                
                if nonempty(1,e) == 1
                    % Now the real deal. For every non-empty epoch, each state type (rem, rippletime,
                    % sleep, sws, and waking) will loop through the times
                    % that that state occurs, counting the number of spikes
                    % that occur during the duration of that occurance. In the end
                    % these spike counts will be summed up and divided by the
                    % total duration of that state to get the mean rate.
                    
                    epochData = C_allStates{s2,r}{1,e}; % Data on a single rat, state, and epoch.
                    nrnsAlltets = [C_allSpikes{1,r}{1,e}{:}]; % Combining nrn spike data from all tets.
                    STs_inState = cell(1,length(nrnsAlltets)); % Cell of spike times that occur within the given state.
                    % Columns are neurons.
                    epochisPyr = zeros(1,size(nrnsAlltets,2));
                    epochisInh = zeros(1,size(nrnsAlltets,2));
        
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
                    % Dur_ratRef(s2,r,e) = r;
                    FRmeans = NaN(size(STs_inState)); % To hold mean firing rate of nrns.
                    for nrn = 1:size(nrnsAlltets,2)
                        if ~isempty(STs_inState{1,nrn}) % Important to leave nrns with no spikes as NaNs.
                            meanFR = length(STs_inState{1,nrn})/stateDur; % Calculate mean firing rate
                            FRmeans(1,nrn) =  meanFR;
                        end

                        % Identify pyramidal vs inhibitory neurons
                        if isfield(nrnsAlltets{nrn},'meanrate')
                            if nrnsAlltets{nrn}.meanrate <= 7 % Call this pyramidal
                                epochisPyr(1,nrn) = 1;
                            elseif nrnsAlltets{nrn}.meanrate > 7 % Call this inhibitory
                                epochisInh(1,nrn) = 1;
                            end
                        end
                    end
                     
                    ratMeans(:,e) = FRmeans; % Stores all FR means for that epoch.
                    occEpochs{1,e} = occTimes;
                    ratisPyr(:,e) = epochisPyr;
                    ratisInh(:,e) = epochisInh;
               
                end
            end
            
            C_stateFRs{s2,a}{1,r} = ratMeans;
            
            FR_allStates{s2,a} = [FR_allStates{s2,a}; ratMeans];    
            Occ_allStates{s2,r} = occEpochs;
            if flag == 0 % This must only be run once, else I will concatenate 
                % the same info onto the matrices for (brain regions) x
                % (states) number of times.

                % ratisPyr and ratisInh mats may have columns that are
                % zeros if the first state doesn't occur during some epochs. This compares
                % the number of epochs a cell is considered either pyramidal or inhibitory
                % and makes a decision based on majority. Ties are discarded.
                for nrn = 1:size(ratisPyr,1)
                    if sum(ratisPyr(nrn,:)) < sum(ratisInh(nrn,:))
                        ratisInh(nrn,:) = 1;
                        ratisPyr(nrn,:) = 0;
                    elseif sum(ratisPyr(nrn,:)) > sum(ratisInh(nrn,:))
                        ratisPyr(nrn,:) = 1;
                        ratisInh(nrn,:) = 0;
                    elseif sum(ratisPyr(nrn,:)) == sum(ratisInh(nrn,:))
                        ratisPyr(nrn,:) = 0;
                        ratisInh(nrn,:) = 0;
                    end
                
                end

                isPyrMat = [isPyrMat; ratisPyr];
                isInhMat = [isInhMat; ratisInh];       
                C_isPyr{1,r} = logical(ratisPyr);
                C_isInh{1,r} = logical(ratisInh);
            end
            if r == size(C_allStates,2) % trips flag after all rats have been run through once.
                flag = 1;
            end
        end
    fprintf("       FRs for state %d completed. \n",s2)
    end
end



% Convert to logical arrays
isPyrMat = logical(isPyrMat);
isInhMat = logical(isInhMat);



end

