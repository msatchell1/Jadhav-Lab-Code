function [C_combStates] = create_combStates(statesMat,stateNames,dataDir)
%CREATE_COMBSTATES Concatenates states and labels them.
%   Combines physiological states for multiple animals in temporal order.
%   Each epoch is handled separately. Also saves a file combinedStates.mat
%   to each animal's folder under the parent directory dataDir.
%
% Inputs:
%
%   statesMat - matrix containing state structures for each state, animal,
%   and epoch. Dimensions (num states) x (num animals), each element
%   containing 1 x (num epochs) arrays, each element containing a struct
%   with fields 'starttime' and 'endtime'. These fields are 1D column
%   vectors with equal length of all the start and end times for that state
%   in that epoch. 
%
%   stateNames - Labels for the different states in statesMat. A
%   column vector with a label for each row of statesMat.
%   
% Outputs:
%
%   C_combStates - a (num animals) x (num epochs) cell array, each
%   element being a struct containing the concatenated and labeled states.


%% Determine which states exist in which epochs. Note this asssumes that all 
% animals have the same state-epoch relationship, as the first animal is
% used as a reference.

% stateExist = cell(size(statesMat{1,1})); % Array to hold the indices of all the existing states 
% % within each epoch.
% 
% for s = 1:size(statesMat,1)
% 
%     for e = 1:size(statesMat{s,1},2)
% 
%         if ~isempty(statesMat{s,1}{1,e})
%             stateExist{1,e} = [stateExist{1,e}; s]; % Appends state to epoch 
%             % if it exists for that epoch.
% 
%         end
%     end
% end


C_combStates = cell(size(statesMat,2),size(statesMat{1,1},2));


% Loop to concatenate the state times
for r = 1:size(statesMat,2)
    
    % Cell matrix (states) x (epochs) with the list of starttimes (col 1)
    % endtimes (col 2) and state index (col 3) within.
    startEndMat = cell(size(statesMat,1),size(statesMat{1,r},2));
    
    for e = 1:size(statesMat{1,r},2)

        % statesData = statesMat(stateExist{1,e},r); % Cell array of relevent states
        
        % Loop through all the states for this epoch getting the start and
        % end times.
        for s = 1:size(statesMat,1)
            if ~isempty(statesMat{s,r}{1,e}) % If state exists for this epoch
                esData = statesMat{s,r}{1,e}; % Epoch state data struct
                if ~isempty(esData)
                    % Assign values times and state ID to matrix
                    startEndMat{s,e} = [esData.starttime,esData.endtime,repmat(s,[size(esData.starttime,1),1])];
                end
            end
        end

        % Combine epoch's state occurances into a single list
        eOccs = vertcat(startEndMat{:,e});
        % Sort based on start times
        [~,sortI] = sort(eOccs(:,1));
        % Reorder list based on sort indices to get a sorted list of the
        % occurances for each state based on when they start.
        sortedOccs = eOccs(sortI,:);

        S_epoch.ratID = r;
        S_epoch.epoch = e;
        S_epoch.stateNames = stateNames';
        S_epoch.sepData = startEndMat(:,e); % State data seperated
        S_epoch.combSortData = sortedOccs; 
        S_epoch.dataDscrp = "Data is sorted into 3 columns. Col 1 = start time (s)," + ...
            " Col 2 = end time, and Col 3 is the index to stateNames.";
        
        C_combStates{r,e} = S_epoch;

    end

    

    


end


end

