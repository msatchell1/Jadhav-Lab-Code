function [C_alldata] = load_data(dataDir,loadRats,filetypes)
%LOAD_DATA Loads .mat data.
%   Detailed explanation goes here



C_alldata = {}; % Cell array to hold data for all rats. If multiple filetypes 
% are loaded, each row holds a different file type, ordered in the same
% order as the elements in filetypes.

disp("Loading new animal data... ")
for r = 1:length(loadRats)
    fprintf("Loading animal: %s \n",loadRats{r})

    for ft = 1:length(filetypes)    

        shortName = loadRats{r};
        chop_idx = strfind(loadRats{r},'_') - 1;
        if ~isempty(chop_idx)
            shortName = loadRats{r}(1:chop_idx); % Gets the first characters of the rat's name before an '_'.
            % So far this is only needed for ER1_NEW to remove the '_NEW'.
        end

        % Does not load from EEG folder
        FileDir = dir(dataDir+"/"+loadRats(r)+'_direct'+"/"+shortName+filetypes{ft}+"*");
    
        if isempty(FileDir)
            error("%s file does not exist for animal: %s \n",filetypes{ft},loadRats{r})
        elseif length(FileDir) > 1
            error("More than one file detected when searching for: %s, in animal: %s \n" + ...
                "Change names in filetypes to be more specific.", filetypes{ft},loadRats{r});
        else
            file = struct2cell(load(string(fullfile(FileDir.folder, FileDir.name)))); % load data
            file = file{:};
            if iscell(file)
                C_alldata{ft,r} = file{1,1};
            else
                C_alldata{ft,r} = file;
            end
            fprintf("       Loaded file: %s  \n",  FileDir.name)
        end
    end
end

C_alldata = clip_17_epochs(C_alldata); % removes extra epoch data.

end

