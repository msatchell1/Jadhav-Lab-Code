% This script will look at raster plots of populations. I want to see SWR
% events in a population of cells. To do this, I will need to sort cells by
% the peaks of their firing fields and plot them in a raster with the SWR
% times marked. 

% Maybe I can also look for theta sequences during the task performance


%% Load data

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.
load_rats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};


filetype = {'linpos01','linfields01'};

C_alldata = {}; % Cell array to hold data for each rat. If multiple filetypes 
% are loaded, each row holds a different file type, ordered in the same
% order as the elements in filetype.

for r = 1:length(load_rats)

    for ft = 1:length(filetype)    

        % Does not load from EEG folder
        File_dir = dir(data_dir+"/"+load_rats(r)+'_direct'+"/*"+filetype{ft}+"*");
    
        if isempty(File_dir)
            error("%s file does not exist for animal: %s \n",filetype{ft},load_rats{r})
        elseif length(File_dir) > 1
            error("More than one file detected when searching for: %s, in animal: %s \n" + ...
                "Must load only one file per element in filetype.", filetype{ft},load_rats{r})
        else
            file = struct2cell(load(string(fullfile(File_dir.folder, File_dir.name)))); % load data
            file = file{:};
            C_alldata{ft,r} = file{1,1};
            fprintf("Loaded animal: %s      file: %s  \n", load_rats{r}, File_dir.name)
        end
    end
end
