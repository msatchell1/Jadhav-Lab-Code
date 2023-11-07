function [] = create_combinedStates(statesMat,stateNames,dataDir)
%CREATE_COMBINEDSTATES Concatenates states and labels them.
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
%   vectors of all the start and end times for that state in that epoch. 
%
%   stateNames - Labels for the different states in statesMat. A
%   column vector with a label for each row of statesMat.
%   
% Outputs:
%
%   S_combinedStates - a (num animals) x (num epochs) cell array, each
%   element being a struct containing the concatenated and labeled states.


%% Determine which states exist in which epochs. Note this asssumes that all 
% animals have the same state-epoch relationship, as the first animal is
% used as a reference.

stateExist = cell(size(statesMat{1,1})); % Array to hold the indices of all the existing states 
% within each epoch.

for s = 1:size(statesMat,1)

    for e = 1:size(statesMat{s,1},2)

        if ~isempty(statesMat{s,1}{1,e})
            stateExist{1,e} = [stateExist{1,e}; s]; % Appends state to epoch 
            % if it exists for that epoch.

        end
    end
end


% Now that we know what states to loop through for each epoch, lets
% concatenate the state times
for r = 1:size(statesMat,2)
    
    for e = 1:size(stateExist,2)

        statesData = statesMat(stateExist{1,e},r); % Cell array of relevent states
        epochData = [statesData{:,1};];

    end
end


end

