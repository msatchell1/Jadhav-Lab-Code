% The same problem that existed for the LFP data, where the timerange needs
% to be shortened to account fo the data that was dropped out during the
% second epoch of ER1. 

%% Load data
% Loads individual tetrode LFP data from aa single or multiple epochs into
% allData

clearvars;

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
load_rats = {'ER1_NEW'};


% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','pos01'

filetypes = 'spikes01'; % only 1 entry allowed

allData = {};

% Only made to work for one filetype
disp("Loading new animal data... ")
for r = 1:length(load_rats)
    fprintf("Loading animal: %s \n",load_rats{r})

    short_name = load_rats{r};
    chop_idx = strfind(load_rats{r},'_') - 1;
    if ~isempty(chop_idx)
        short_name = load_rats{r}(1:chop_idx); % Gets the first characters of the rat's name before an '_'.
        % So far this is only needed for ER1_NEW to remove the '_NEW'.
    end


    File_dir = dir(data_dir+"/"+load_rats(r)+'_direct'+"/"+short_name+filetypes+"*");

    if isempty(File_dir)
        error("%s file does not exist for animal: %s \n",filetypes,load_rats{r})
    % elseif length(File_dir) > 1
    %     error("More than one file detected when searching for: %s, in animal: %s \n" + ...
    %         "Change names in filetypes to be more specific.", filetypes{ft},load_rats{r});
    else
        for f = 1:length(File_dir)
            file = struct2cell(load(string(fullfile(File_dir(f).folder, File_dir(f).name)))); % load data
            file = file{:};
            allData{f,r} = file{1,1};
            fprintf("       Loaded file: %s  \n",  File_dir.name)
        end
    end
end


%% Edit timerange values

for f = 1:size(allData,1)
    for e = 1:size(allData{f,1},2)
        tets = allData{f,1}{1,e};
        for t = 1:size(tets,2)
            if ~isempty(tets{1,t})
                for n = 1:size(tets(t),2) % Loops through all neurons on a tetrode
                    if isfield(tets{1,t}{1,n}, 'timerange')
                        S = tets{1,t}{1,n};
                        timeVals = S.timerange(1) : 1/1500 : S.timerange(2);
                        dpDiff = 460; % data point difference, i.e. discrepancy
                        % between the amount of data that should be there according to timerange
                        % and the actual amount of data.
                        timeDiff = dpDiff/1500; 
                        % Make the timerange shorter
                        allData{f,1}{1,e}{1,t}{1,n}.timerange(2) = S.timerange(2) - timeDiff;
                    end
                end
            end
        end
    end
end



%% Save the edited structs back to their original location

for f = 1:length(File_dir)
    % This is SO stupid, but because matlab's load/save functions retain the
    % name of the struct thsat was saved, in order for my code to work
    % the struct names must match those used by Justin, so I need a
    % different struct name for each of the EEG filetypes.
    if contains(filetypes,"spikes")
        spikes = allData(f,1);
        save(string(fullfile(File_dir(f).folder, File_dir(f).name)), "spikes");
    end
end