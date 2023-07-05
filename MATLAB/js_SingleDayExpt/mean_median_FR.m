% This script will recreate the results that Aanchal got on median and mean
% firing rates of neurons in CA1 and PFC during awake and sleep.
% Michael Satchell 07/03/2023


%% Load data

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.
load_rats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Mean firing rate data for each nrn is held in the cellinfo file.
filetype = {'cellinfo'};

C_alldata = {}; % Cell array to hold data for each rat

for i = 1:length(load_rats)

    File_dir = dir(data_dir+"/"+load_rats(i)+'_direct'+"/*"+filetype(1)+"*");
    % There should only be one cellinfo file per animal
    if isempty(File_dir)
        error("cellinfo file does not exist for animal: %s \n",load_rats{i})
    elseif length(File_dir) > 1
        error("Only 1 file per animal expected. Number loaded for animal %s: %d \n", load_rats{i},lenth(File_dir))
    else
        file = struct2cell(load(string(fullfile(File_dir.folder, File_dir.name)))); % load data
        file = file{:};
        C_alldata{i} = file{1,1};
        fprintf("Loaded animal: %s      Num epochs: %d      Num cells: ~%d \n", load_rats{i},length(file{1,1}),length([file{1,1}{1,1}{:}]))
        % Note that the number of cells on a tetrode varies across the
        % epochs, so the number listed here is only an estimate based on
        % the first epoch.
    end
end



%% Plot mean and median of mean firing rates for one rat

animalnum = 2;
ratdata = C_alldata{1,animalnum}; % this works for any animal by just indexing the animalnum into C_alldata.

% This is to get the largest number of neurons in an epoch. row 1 for CA1
% nrns, row 2 for PFC.
nrn_nums = zeros(1,length(ratdata));
for e = 1:length(ratdata) % Loop through epochs
    flat_nrns = [ratdata{1,e}{:}];
    for n = 1:length(flat_nrns)
        nrn_nums(e) = nrn_nums(e) + 1;
    end
end
CA1_mat = NaN([max(nrn_nums),length(ratdata)]); % Firing rate matrix. I CANNOT assume
% that a given row of this matrix holds the same cell across epochs,
% because the number of cells on tetrodes changes every epoch.
PFC_mat = NaN([max(nrn_nums),length(ratdata)]);

% Now populate matrix with mean firing rates for every nrn in every epoch.
for e = 1:length(ratdata)

    flat_nrns = [ratdata{1,e}{:}];
    for n = 1:length(flat_nrns)
        if isfield(flat_nrns{1,n}, 'meanrate')
            if strcmp(flat_nrns{1,n}.area, 'CA1')
                CA1_mat(n,e) = flat_nrns{1,n}.meanrate;
            elseif strcmp(flat_nrns{1,n}.area, 'PFC')
                PFC_mat(n,e) = flat_nrns{1,n}.meanrate;
            end
        end
    end
end



figure;
sgtitle(sprintf("Animal %s, Epoch 1 CA1 Cells: %d, Epoch 1 PFC Cells: %d",...
    load_rats{animalnum},sum(~isnan(CA1_mat(:,1))), sum(~isnan(PFC_mat(:,1)))), 'Interpreter','none')
% I need to read up on statistics to figure out if I should be plotting the
% standard deviation of the mean firing rates for the error bars, or
% something else.
subplot(2,2,1)
title("CA1 Mean FR")
ylabel("Mean Firing Rate (Hz)")
xlabel("Epoch")
hold on;
sleep_inds = 1:2:length(ratdata);
wake_inds = 2:2:length(ratdata);
CA1_std = std(CA1_mat,0,1,"omitnan");
CA1_mean = mean(CA1_mat,1,"omitnan");
CA1_sleep_std = CA1_std(sleep_inds);
CA1_wake_std = CA1_std(wake_inds);
CA1_sleep_mean = CA1_mean(sleep_inds);
CA1_wake_mean = CA1_mean(wake_inds);
errorbar(sleep_inds, CA1_sleep_mean, CA1_sleep_std, 'o',Color=[0 0.4470 0.7410], MarkerFaceColor='blue')
errorbar(wake_inds, CA1_wake_mean, CA1_wake_std, 'o', Color=[0 0.4470 0.7410])
legend("sleep","wake");

subplot(2,2,2)
title("PFC Mean FR")
ylabel("Mean Firing Rate (Hz)")
xlabel("Epoch")
hold on;
PFC_std = std(PFC_mat,0,1,"omitnan");
PFC_mean = mean(PFC_mat,1,"omitnan");
PFC_sleep_std = PFC_std(sleep_inds);
PFC_wake_std = PFC_std(wake_inds);
PFC_sleep_mean = PFC_mean(sleep_inds);
PFC_wake_mean = PFC_mean(wake_inds);
errorbar(sleep_inds, PFC_sleep_mean, PFC_sleep_std, 'o',Color=[0.8500 0.3250 0.0980], MarkerFaceColor=[0.8500 0.3250 0.0580])
errorbar(wake_inds, PFC_wake_mean, PFC_wake_std, 'o', Color=[0.8500 0.3250 0.0980])
legend("sleep","wake");

subplot(2,2,3)
title("CA1 Median FR")
ylabel("Median Firing Rate (Hz)")
xlabel("Epoch")
hold on;
CA1_std = std(CA1_mat,0,1,"omitnan");
CA1_median = median(CA1_mat,1,"omitnan");
CA1_sleep_std = CA1_std(sleep_inds);
CA1_wake_std = CA1_std(wake_inds);
CA1_sleep_median = CA1_median(sleep_inds);
CA1_wake_median = CA1_median(wake_inds);
errorbar(sleep_inds, CA1_sleep_median, CA1_sleep_std, 'o',Color=[0 0.4470 0.7410], MarkerFaceColor='blue')
errorbar(wake_inds, CA1_wake_median, CA1_wake_std, 'o', Color=[0 0.4470 0.7410])
legend("sleep","wake");

subplot(2,2,4)
title("PFC Median FR")
ylabel("Median Firing Rate (Hz)")
xlabel("Epoch")
hold on;
PFC_std = std(PFC_mat,0,1,"omitnan");
PFC_median = median(PFC_mat,1,"omitnan");
PFC_sleep_std = PFC_std(sleep_inds);
PFC_wake_std = PFC_std(wake_inds);
PFC_sleep_median = PFC_median(sleep_inds);
PFC_wake_median = PFC_median(wake_inds);
errorbar(sleep_inds, PFC_sleep_median, PFC_sleep_std, 'o',Color=[0.8500 0.3250 0.0980], MarkerFaceColor=[0.8500 0.3250 0.0580])
errorbar(wake_inds, PFC_wake_median, PFC_wake_std, 'o', Color=[0.8500 0.3250 0.0980])
legend("sleep","wake");



%% Plotting Across Animals

% Figure out animal with longest number of epochs (should be 19)
epochcounts = zeros([1,length(C_alldata)]);
for a = 1:length(C_alldata)
    epochcounts(a) = length(C_alldata{1,a});
end

CA1_std_mat = NaN(max(epochcounts),length(C_alldata));
CA1_mean_mat = NaN(max(epochcounts),length(C_alldata));
CA1_median_mat = NaN(max(epochcounts),length(C_alldata));
PFC_std_mat = NaN(max(epochcounts),length(C_alldata));
PFC_mean_mat = NaN(max(epochcounts),length(C_alldata));
PFC_median_mat = NaN(max(epochcounts),length(C_alldata));

% Loop all animals
for a = 1:length(C_alldata)


    ratdata = C_alldata{1,a};
    
    % This is to get the largest number of neurons in an epoch. row 1 for CA1
    % nrns, row 2 for PFC.
    nrn_nums = zeros(1,length(ratdata));
    for e = 1:length(ratdata) % Loop through epochs
        flat_nrns = [ratdata{1,e}{:}];
        for n = 1:length(flat_nrns)
            nrn_nums(e) = nrn_nums(e) + 1;
        end
    end
    CA1_mat = NaN([max(nrn_nums),length(ratdata)]); % Firing rate matrix. I CANNOT assume
    % that a given row of this matrix holds the same cell across epochs,
    % because the number of cells on tetrodes changes every epoch.
    PFC_mat = NaN([max(nrn_nums),length(ratdata)]);
    
    % Now populate matrix with mean firing rates for every nrn in every epoch.
    for e = 1:length(ratdata)
    
        flat_nrns = [ratdata{1,e}{:}];
        for n = 1:length(flat_nrns)
            if isfield(flat_nrns{1,n}, 'meanrate')
                if strcmp(flat_nrns{1,n}.area, 'CA1')
                    CA1_mat(n,e) = flat_nrns{1,n}.meanrate;
                elseif strcmp(flat_nrns{1,n}.area, 'PFC')
                    PFC_mat(n,e) = flat_nrns{1,n}.meanrate;
                end
            end
        end
    end


    CA1_std = std(CA1_mat,0,1,"omitnan");
    CA1_mean = mean(CA1_mat,1,"omitnan");
    CA1_median = median(CA1_mat,1,"omitnan");
    PFC_std = std(PFC_mat,0,1,"omitnan");
    PFC_mean = mean(PFC_mat,1,"omitnan");
    PFC_median = median(PFC_mat,1,"omitnan");

    % Now to add the data to the matrices, I need to add NaN elements to
    % the vectors so they can fill in the matrix columns.
    if length(ratdata) < max(epochcounts)
        CA1_std(end+1:max(epochcounts)) = nan;
        CA1_mean(end+1:max(epochcounts)) = nan;
        CA1_median(end+1:max(epochcounts)) = nan;
        PFC_std(end+1:max(epochcounts)) = nan;
        PFC_mean(end+1:max(epochcounts)) = nan;
        PFC_median(end+1:max(epochcounts)) = nan;
    end


    CA1_std_mat(:,a) = CA1_std';
    CA1_mean_mat(:,a) = CA1_mean';
    CA1_median_mat(:,a) = CA1_median';
    PFC_std_mat(:,a) = PFC_std';
    PFC_mean_mat(:,a) = PFC_mean';
    PFC_median_mat(:,a) = PFC_median';

end


% I'm not sure whether to use the mean or median data for the aggregate
% analyses, so Im going to go with mean for now.
% Here the average across all animals is taken to get a mean firing rate
% across all epochs.

a_CA1_mean = mean(CA1_mean_mat,2,'omitnan'); % mean across animals
a_PFC_mean = mean(PFC_mean_mat,2,'omitnan'); 

a_CA1_std = std(CA1_mean_mat,0,2,'omitnan'); % Std across animals
a_PFC_std = std(PFC_mean_mat,0,2,'omitnan');

% Indices for sleep and wake
sleep_inds = 1:2:size(a_CA1_mean,1);
wake_inds = 2:2:size(a_CA1_mean,1);

figure;
sgtitle(sprintf("Mean Firing Rates Across %d Animals",length(load_rats)));

subplot(1,2,1)
hold on;
title("CA1")
ylabel("Mean Firing Rate (Hz)")
xlabel("Epoch")
errorbar(sleep_inds, a_CA1_mean(sleep_inds), a_CA1_std(sleep_inds), 'o',Color=[0 0.4470 0.7410], MarkerFaceColor='blue')
errorbar(wake_inds, a_CA1_mean(wake_inds), a_CA1_std(wake_inds), 'o', Color=[0 0.4470 0.7410])
legend("sleep","wake");

subplot(1,2,2)
hold on;
title("PFC")
ylabel("Mean Firing Rate (Hz)")
xlabel("Epoch")
errorbar(sleep_inds, a_PFC_mean(sleep_inds), a_PFC_std(sleep_inds), 'o',Color=[0.8500 0.3250 0.0980], MarkerFaceColor=[0.8500 0.3250 0.0580])
errorbar(wake_inds, a_PFC_mean(wake_inds), a_PFC_std(wake_inds), 'o', Color=[0.8500 0.3250 0.0980])
legend("sleep","wake");





