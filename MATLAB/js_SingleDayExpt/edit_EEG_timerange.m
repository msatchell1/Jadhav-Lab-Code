% Edit the timerange field on the eeg file structs so that they match the
% actual length of data in the data field.


%% Load data from EEG folder
% Loads individual tetrode LFP data from aa single or multiple epochs into
% allData

clearvars;

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
load_rats = {'ER1_NEW'};

eStr = "02"; % epoch number
e = str2num(eStr);


% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','pos01'

filetypes = 'ripple01-'+eStr; % only 1 entry allowed

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


    File_dir = dir(data_dir+"/"+load_rats(r)+'_direct'+"/EEG/"+short_name+filetypes+"*");

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

for t = 1:size(allData,1) % Note t might not be the tetrode number if any tetrodes are missing from the EEG folder
    tets = allData{t,1}{1,e};
    tetnum = find(~cellfun(@isempty,tets)); % actual tetrode number
    S = tets{1,tetnum};
    timeVals = S.timerange(1) : 1/S.samprate : S.timerange(2);
    dpDiff = length(timeVals) - length(S.data); % data point difference, i.e. discrepancy
    % between the amount of data that should be there according to timerange
    % and the actual amount of data.
    timeDiff = dpDiff/S.samprate; 
    % Make the timerange shorter
    allData{t,1}{1,e}{1,tetnum}.timerange(2) = S.timerange(2) - timeDiff;
    % Also important to update 'endtime' field
    allData{t,1}{1,e}{1,tetnum}.endtime = S.timerange(2) - timeDiff;
end




%% Save the edited structs back to their original location

for f = 1:length(File_dir)
    % This is SO stupid, but because matlab's load/save functions retain the
    % name of the struct thsat was saved, in order for my code to work
    % the struct names must match those used by Justin, so I need a
    % different struct name for each of the EEG filetypes.
    if contains(filetypes,"eegref01")
        eegref = allData(f,1);
        save(string(fullfile(File_dir(f).folder, File_dir(f).name)), "eegref");
    elseif contains(filetypes,"eeg01")
        eeg = allData(f,1);
        save(string(fullfile(File_dir(f).folder, File_dir(f).name)), "eeg");
    elseif contains(filetypes,"delta01")
        delta = allData(f,1);
        save(string(fullfile(File_dir(f).folder, File_dir(f).name)), "delta");
    elseif contains(filetypes,"gamma01")
        gamma = allData(f,1);
        save(string(fullfile(File_dir(f).folder, File_dir(f).name)), "gamma");
    elseif contains(filetypes,"ripple01")
        ripple = allData(f,1);
        save(string(fullfile(File_dir(f).folder, File_dir(f).name)), "ripple");
    elseif contains(filetypes,"theta01")
        theta = allData(f,1);
        save(string(fullfile(File_dir(f).folder, File_dir(f).name)), "theta");
    end
end