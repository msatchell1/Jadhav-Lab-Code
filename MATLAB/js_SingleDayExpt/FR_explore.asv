% This script will recreate some results that Aanchal got on looking at
% firing rates of neurons in CA1 and PFC during awake and sleep. It will
% also explore firing rate statistics in general.
% Michael Satchell 07/03/2023


%% Load data

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
load_rats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common filetypes: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01'.
filetype = {'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01'};

C_alldata = {}; % Cell array to hold data for all rats. If multiple filetypes 
% are loaded, each row holds a different file type, ordered in the same
% order as the elements in filetype.

disp("Loading new animal data... ")
for r = 1:length(load_rats)
    fprintf("Loaded animal: %s \n",load_rats{r})

    for ft = 1:length(filetype)    

        short_name = load_rats{r};
        chop_idx = strfind(load_rats{r},'_') - 1;
        if ~isempty(chop_idx)
            short_name = load_rats{r}(1:chop_idx); % Gets the first characters of the rat's name before an '_'.
            % So far this is only needed for ER1_NEW to remove the '_NEW'.
        end

        % Does not load from EEG folder
        File_dir = dir(data_dir+"/"+load_rats(r)+'_direct'+"/"+short_name+filetype{ft}+"*");
    
        if isempty(File_dir)
            error("%s file does not exist for animal: %s \n",filetype{ft},load_rats{r})
        elseif length(File_dir) > 1
            error("More than one file detected when searching for: %s, in animal: %s \n" + ...
                "Change names in filetype to be more specific.", filetype{ft},load_rats{r});
        else
            file = struct2cell(load(string(fullfile(File_dir.folder, File_dir.name)))); % load data
            file = file{:};
            C_alldata{ft,r} = file{1,1};
            fprintf("       Loaded file: %s  \n",  File_dir.name)
        end
    end
end


% load_rats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};
% 
% % Mean firing rate data for each nrn is held in the cellinfo file.
% filetype = {'cellinfo'};
% 
% C_alldata = {}; % Cell array to hold data for each rat
% 
% for i = 1:length(load_rats)
% 
%     File_dir = dir(data_dir+"/"+load_rats(i)+'_direct'+"/*"+filetype(1)+"*");
%     % There should only be one cellinfo file per animal
%     if isempty(File_dir)
%         error("cellinfo file does not exist for animal: %s \n",load_rats{i})
%     elseif length(File_dir) > 1
%         error("Only 1 file per animal expected. Number loaded for animal %s: %d \n", load_rats{i},lenth(File_dir))
%     else
%         file = struct2cell(load(string(fullfile(File_dir.folder, File_dir.name)))); % load data
%         file = file{:};
%         C_alldata{i} = file{1,1};
%         fprintf("Loaded animal: %s      Num epochs: %d      Num cells: ~%d \n", load_rats{i},length(file{1,1}),length([file{1,1}{1,1}{:}]))
%         % Note that the number of cells on a tetrode varies across the
%         % epochs, so the number listed here is only an estimate based on
%         % the first epoch.
%     end
% end



%% Plot mean and median of mean firing rates for one rat

% animalnum = 2;
% ratdata = C_alldata{1,animalnum}; % this works for any animal by just indexing the animalnum into C_alldata.
% 
% % This is to get the largest number of neurons in an epoch. row 1 for CA1
% % nrns, row 2 for PFC.
% nrn_nums = zeros(1,length(ratdata));
% for e = 1:length(ratdata) % Loop through epochs
%     flat_nrns = [ratdata{1,e}{:}];
%     for n = 1:length(flat_nrns)
%         nrn_nums(e) = nrn_nums(e) + 1;
%     end
% end
% CA1_mat = NaN([max(nrn_nums),length(ratdata)]); % Firing rate matrix. I CANNOT assume
% % that a given row of this matrix holds the same cell across epochs,
% % because the number of cells on tetrodes changes every epoch.
% PFC_mat = NaN([max(nrn_nums),length(ratdata)]);
% CA1_numnrns = zeros(1,length(ratdata));
% PFC_numnrns = zeros(1,length(ratdata));
% 
% % Now populate matrix with mean firing rates for every nrn in every epoch.
% for e = 1:length(ratdata)
% 
%     flat_nrns = [ratdata{1,e}{:}];
%     for n = 1:length(flat_nrns)
%         if isfield(flat_nrns{1,n}, 'meanrate')
%             if strcmp(flat_nrns{1,n}.area, 'CA1')
%                 CA1_mat(n,e) = flat_nrns{1,n}.meanrate;
%                 CA1_numnrns(e) = CA1_numnrns(e) + 1;
%             elseif strcmp(flat_nrns{1,n}.area, 'PFC')
%                 PFC_mat(n,e) = flat_nrns{1,n}.meanrate;
%                 PFC_numnrns(e) = PFC_numnrns(e) + 1;
%             end
%         end
%     end
% end
% 
% % Later I want to calculate confidence intervals for plotting error
% CA1_CI95 = {}; % To store probability
% PFC_CI95 = {};
% for e = 1:size(CA1_mat,2)
%     cCI95 = tinv([0.025 0.975], CA1_numnrns(e)-1);
%     pCI95 = tinv([0.025 0.975], PFC_numnrns(e)-1);
%     CA1_CI95{e} = cCI95;
%     PFC_CI95{e} = pCI95;
% end
% 
% figure;
% sgtitle(sprintf("Animal %s, Epoch 1 CA1 Cells: %d, Epoch 1 PFC Cells: %d",...
%     load_rats{animalnum},sum(~isnan(CA1_mat(:,1))), sum(~isnan(PFC_mat(:,1)))), 'Interpreter','none')
% % I need to read up on statistics to figure out if I should be plotting the
% % standard deviation of the mean firing rates for the error bars, or
% % something else.
% subplot(2,2,1)
% title("CA1 Mean FR (95% CI Error)")
% ylabel("Mean Firing Rate (Hz)")
% xlabel("Epoch")
% hold on;
% sleep_inds = 1:2:length(ratdata);
% wake_inds = 2:2:length(ratdata);
% CA1_std = std(CA1_mat,0,1,"omitnan");
% CA1_mean = mean(CA1_mat,1,"omitnan");
% CA1_sleep_std = CA1_std(sleep_inds);
% CA1_wake_std = CA1_std(wake_inds);
% CA1_sleep_mean = CA1_mean(sleep_inds);
% CA1_wake_mean = CA1_mean(wake_inds);
% 
% CA1_yCI95 = {}; % Calculating confidence interval error
% for e = 1:length(ratdata)
%     ySEM = CA1_std./sqrt(CA1_numnrns);
%     CA1_yCI95{e} = bsxfun(@times, ySEM(e), CA1_CI95{1,e}(:)); 
% end
% CA1_yCI95 = [CA1_yCI95{:}];
% 
% errorbar(sleep_inds, CA1_sleep_mean, CA1_yCI95(1,sleep_inds), CA1_yCI95(2,sleep_inds), 'o',Color=[0 0.4470 0.7410], MarkerFaceColor='blue')
% errorbar(wake_inds, CA1_wake_mean, CA1_yCI95(1,wake_inds), CA1_yCI95(2,wake_inds), 'o', Color=[0 0.4470 0.7410])
% legend("sleep","wake");
% 
% subplot(2,2,2)
% title("PFC Mean FR")
% ylabel("Mean Firing Rate (Hz)")
% xlabel("Epoch")
% hold on;
% PFC_std = std(PFC_mat,0,1,"omitnan");
% PFC_mean = mean(PFC_mat,1,"omitnan");
% PFC_sleep_std = PFC_std(sleep_inds);
% PFC_wake_std = PFC_std(wake_inds);
% PFC_sleep_mean = PFC_mean(sleep_inds);
% PFC_wake_mean = PFC_mean(wake_inds);
% errorbar(sleep_inds, PFC_sleep_mean, PFC_sleep_std, 'o',Color=[0.8500 0.3250 0.0980], MarkerFaceColor=[0.8500 0.3250 0.0580])
% errorbar(wake_inds, PFC_wake_mean, PFC_wake_std, 'o', Color=[0.8500 0.3250 0.0980])
% legend("sleep","wake");
% 
% subplot(2,2,3)
% title("CA1 Median FR")
% ylabel("Median Firing Rate (Hz)")
% xlabel("Epoch")
% hold on;
% CA1_std = std(CA1_mat,0,1,"omitnan");
% CA1_median = median(CA1_mat,1,"omitnan");
% CA1_sleep_std = CA1_std(sleep_inds);
% CA1_wake_std = CA1_std(wake_inds);
% CA1_sleep_median = CA1_median(sleep_inds);
% CA1_wake_median = CA1_median(wake_inds);
% errorbar(sleep_inds, CA1_sleep_median, CA1_sleep_std, 'o',Color=[0 0.4470 0.7410], MarkerFaceColor='blue')
% errorbar(wake_inds, CA1_wake_median, CA1_wake_std, 'o', Color=[0 0.4470 0.7410])
% legend("sleep","wake");
% 
% subplot(2,2,4)
% title("PFC Median FR")
% ylabel("Median Firing Rate (Hz)")
% xlabel("Epoch")
% hold on;
% PFC_std = std(PFC_mat,0,1,"omitnan");
% PFC_median = median(PFC_mat,1,"omitnan");
% PFC_sleep_std = PFC_std(sleep_inds);
% PFC_wake_std = PFC_std(wake_inds);
% PFC_sleep_median = PFC_median(sleep_inds);
% PFC_wake_median = PFC_median(wake_inds);
% errorbar(sleep_inds, PFC_sleep_median, PFC_sleep_std, 'o',Color=[0.8500 0.3250 0.0980], MarkerFaceColor=[0.8500 0.3250 0.0580])
% errorbar(wake_inds, PFC_wake_median, PFC_wake_std, 'o', Color=[0.8500 0.3250 0.0980])
% legend("sleep","wake");



%% Plotting Across Animals

% % Figure out animal with longest number of epochs (should be 19)
% epochcounts = zeros([1,length(C_alldata)]);
% for a = 1:length(C_alldata)
%     epochcounts(a) = length(C_alldata{1,a});
% end
% 
% CA1_std_mat = NaN(max(epochcounts),length(C_alldata));
% CA1_mean_mat = NaN(max(epochcounts),length(C_alldata));
% CA1_median_mat = NaN(max(epochcounts),length(C_alldata));
% PFC_std_mat = NaN(max(epochcounts),length(C_alldata));
% PFC_mean_mat = NaN(max(epochcounts),length(C_alldata));
% PFC_median_mat = NaN(max(epochcounts),length(C_alldata));
% 
% CA1_allrates = {}; PFC_allrates = {};
% 
% % Loop all animals
% for a = 1:length(C_alldata)
% 
% 
%     ratdata = C_alldata{1,a};
% 
%     % This is to get the largest number of neurons in an epoch. row 1 for CA1
%     % nrns, row 2 for PFC.
%     nrn_nums = zeros(1,length(ratdata));
%     for e = 1:length(ratdata) % Loop through epochs
%         flat_nrns = [ratdata{1,e}{:}];
%         for n = 1:length(flat_nrns)
%             nrn_nums(e) = nrn_nums(e) + 1;
%         end
%     end
%     CA1_mat = NaN([max(nrn_nums),length(ratdata)]); % Firing rate matrix. I CANNOT assume
%     % that a given row of this matrix holds the same cell across epochs,
%     % because the number of cells on tetrodes changes every epoch.
%     PFC_mat = NaN([max(nrn_nums),length(ratdata)]);
% 
%     % Now populate matrix with mean firing rates for every nrn in every epoch.
%     for e = 1:length(ratdata)
% 
%         flat_nrns = [ratdata{1,e}{:}];
%         for n = 1:length(flat_nrns)
%             if isfield(flat_nrns{1,n}, 'meanrate')
%                 if strcmp(flat_nrns{1,n}.area, 'CA1')
%                     CA1_mat(n,e) = flat_nrns{1,n}.meanrate;
%                 elseif strcmp(flat_nrns{1,n}.area, 'PFC')
%                     PFC_mat(n,e) = flat_nrns{1,n}.meanrate;
%                 end
%             end
%         end
%     end
% 
%     % Store all mean firing rates in cell arrays for later use
%     CA1_allrates(a) = {CA1_mat}; PFC_allrates(a) = {PFC_mat};
% 
%     CA1_std = std(CA1_mat,0,1,"omitnan");
%     CA1_mean = mean(CA1_mat,1,"omitnan");
%     CA1_median = median(CA1_mat,1,"omitnan");
%     PFC_std = std(PFC_mat,0,1,"omitnan");
%     PFC_mean = mean(PFC_mat,1,"omitnan");
%     PFC_median = median(PFC_mat,1,"omitnan");
% 
%     % Now to add the data to the matrices, I need to add NaN elements to
%     % the vectors so they can fill in the matrix columns.
%     if length(ratdata) < max(epochcounts)
%         CA1_std(end+1:max(epochcounts)) = nan;
%         CA1_mean(end+1:max(epochcounts)) = nan;
%         CA1_median(end+1:max(epochcounts)) = nan;
%         PFC_std(end+1:max(epochcounts)) = nan;
%         PFC_mean(end+1:max(epochcounts)) = nan;
%         PFC_median(end+1:max(epochcounts)) = nan;
%     end
% 
% 
%     CA1_std_mat(:,a) = CA1_std';
%     CA1_mean_mat(:,a) = CA1_mean';
%     CA1_median_mat(:,a) = CA1_median';
%     PFC_std_mat(:,a) = PFC_std';
%     PFC_mean_mat(:,a) = PFC_mean';
%     PFC_median_mat(:,a) = PFC_median';
% 
% end
% 
% 
% % I'm not sure whether to use the mean or median data for the aggregate
% % analyses, so Im going to go with mean for now.
% % Here the average across all animals is taken to get a mean firing rate
% % across all epochs.
% 
% a_CA1_mean = mean(CA1_mean_mat,2,'omitnan'); % mean across animals
% a_PFC_mean = mean(PFC_mean_mat,2,'omitnan'); 
% a_CA1_median = median(CA1_mean_mat,2,'omitnan'); %median
% a_PFC_median = median(PFC_mean_mat,2,'omitnan');
% 
% a_CA1_std = std(CA1_mean_mat,0,2,'omitnan'); % Std across animals
% a_PFC_std = std(PFC_mean_mat,0,2,'omitnan');
% 
% % Indices for sleep and wake
% sleep_inds = 1:2:size(a_CA1_mean,1);
% wake_inds = 2:2:size(a_CA1_mean,1);

% figure; %mean
% sgtitle(sprintf("Mean Firing Rates Across %d Animals",length(load_rats)));

% subplot(1,2,1)
% hold on;
% title("CA1")
% ylabel("Mean Firing Rate (Hz)")
% xlabel("Epoch")
% errorbar(sleep_inds, a_CA1_mean(sleep_inds), a_CA1_std(sleep_inds), 'o',Color=[0 0.4470 0.7410], MarkerFaceColor='blue')
% errorbar(wake_inds, a_CA1_mean(wake_inds), a_CA1_std(wake_inds), 'o', Color=[0 0.4470 0.7410])
% legend("sleep","wake");
% 
% subplot(1,2,2)
% hold on;
% title("PFC")
% ylabel("Mean Firing Rate (Hz)")
% xlabel("Epoch")
% errorbar(sleep_inds, a_PFC_mean(sleep_inds), a_PFC_std(sleep_inds), 'o',Color=[0.8500 0.3250 0.0980], MarkerFaceColor=[0.8500 0.3250 0.0580])
% errorbar(wake_inds, a_PFC_mean(wake_inds), a_PFC_std(wake_inds), 'o', Color=[0.8500 0.3250 0.0980])
% legend("sleep","wake");
% 
% 
% figure; %median
% sgtitle(sprintf("Median Firing Rates Across %d Animals",length(load_rats)));
% 
% subplot(1,2,1)
% hold on;
% title("CA1")
% ylabel("Median Firing Rate (Hz)")
% xlabel("Epoch")
% errorbar(sleep_inds, a_CA1_median(sleep_inds), a_CA1_std(sleep_inds), 'o',Color=[0 0.4470 0.7410], MarkerFaceColor='blue')
% errorbar(wake_inds, a_CA1_median(wake_inds), a_CA1_std(wake_inds), 'o', Color=[0 0.4470 0.7410])
% legend("sleep","wake");
% 
% subplot(1,2,2)
% hold on;
% title("PFC")
% ylabel("Median Firing Rate (Hz)")
% xlabel("Epoch")
% errorbar(sleep_inds, a_PFC_median(sleep_inds), a_PFC_std(sleep_inds), 'o',Color=[0.8500 0.3250 0.0980], MarkerFaceColor=[0.8500 0.3250 0.0580])
% errorbar(wake_inds, a_PFC_median(wake_inds), a_PFC_std(wake_inds), 'o', Color=[0.8500 0.3250 0.0980])
% legend("sleep","wake");
% 


%% Trying the same thing but treating as if all nrns came from 1 animal. Also plot histograms

% Get a list of the number of epochs across animals.
epochcounts = zeros([1,length(C_alldata)]);
for a = 1:length(C_alldata)
    epochcounts(a) = length(C_alldata{1,a});
end

% First I need to assemble the firing rate information from all the rats
% into one matrix for CA1 and one for PFC.
CA1_allrates = {}; PFC_allrates = {};

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
    
    % Store all mean firing rates in cell arrays for later use
    CA1_allrates(a) = {CA1_mat}; PFC_allrates(a) = {PFC_mat};

end


CA1_1a = []; PFC_1a = []; % "One animal" matrices. Each column will be an epoch,
% and the rows will be appended mean rate data from all rats.

for r = 1:size(CA1_allrates,2)

    MCA1 = CA1_allrates{1,r};
    MPFC = PFC_allrates{1,r};
    if size(MCA1,2) < max(epochcounts)
        % Horizontally append the necessary number of columns of NaNs to
        % the matrix for later appending.
        MCA1 = [MCA1 NaN(size(MCA1,1),(max(epochcounts)-size(MCA1,2)))];
        MPFC = [MPFC NaN(size(MPFC,1),(max(epochcounts)-size(MPFC,2)))];
    end

    CA1_1a = [CA1_1a; MCA1]; PFC_1a = [PFC_1a; MPFC];
end




CA1_numnrns = sum(~isnan(CA1_1a),1);
PFC_numnrns = sum(~isnan(PFC_1a),1);

CA1_1a_mean = mean(CA1_1a,1,'omitnan'); % mean across animals
PFC_1a_mean = mean(PFC_1a,1,'omitnan'); 

CA1_1a_std = std(CA1_1a,0,1,'omitnan'); % Std across animals
PFC_1a_std = std(PFC_1a,0,1,'omitnan');

CA1_1a_SEM = CA1_1a_std./sqrt(CA1_numnrns);
PFC_1a_SEM = PFC_1a_std./sqrt(PFC_numnrns);
CA1_1a_CI = zeros(2,size(CA1_1a_mean,2)); % Note size(CA1_1a_mean,2) is the number of epochs
PFC_1a_CI = zeros(2,size(CA1_1a_mean,2));

CI_val = 95; % Confidence interval to be used in %.

for e = 1:size(CA1_1a_mean,2)
    % CI from a normal? distribution with N-1 data points
    temp_CI_CA1 = tinv([(1-CI_val/100)/2, 1-(1-CI_val/100)/2], CA1_numnrns(1,e)-1); 
    temp_CI_PFC = tinv([(1-CI_val/100)/2, 1-(1-CI_val/100)/2], PFC_numnrns(1,e)-1);

    % multiplying this CI by the SEM to get the error size
    CA1_1a_CI(:,e) = bsxfun(@times, CA1_1a_SEM(e), temp_CI_CA1(:));
    PFC_1a_CI(:,e) = bsxfun(@times, PFC_1a_SEM(e), temp_CI_PFC(:));
end

% Indices for sleep and wake
sleep_inds = 1:2:size(CA1_1a_mean,2);
wake_inds = 2:2:size(CA1_1a_mean,2);

figure;
sgtitle(sprintf("1a Mean Firing Rates Across %d Animals",length(load_rats)));

subplot(1,2,1)
hold on;
title("CA1")
ylabel("Mean Firing Rate (Hz)")
xlabel("Epoch")
errorbar(sleep_inds, CA1_1a_mean(sleep_inds), CA1_1a_CI(1,sleep_inds),CA1_1a_CI(2,sleep_inds), ...
    'o',Color=[0 0.4470 0.7410], MarkerFaceColor='blue')
errorbar(wake_inds, CA1_1a_mean(wake_inds), CA1_1a_CI(1,wake_inds), CA1_1a_CI(2,wake_inds), ...
    'o', Color=[0 0.4470 0.7410])
legend("sleep","wake");

subplot(1,2,2)
hold on;
title("PFC")
ylabel("Mean Firing Rate (Hz)")
xlabel("Epoch")
errorbar(sleep_inds, PFC_1a_mean(sleep_inds), PFC_1a_CI(1,sleep_inds),PFC_1a_CI(2,sleep_inds), ...
    'o',Color=[0.8500 0.3250 0.0980], MarkerFaceColor=[0.8500 0.3250 0.0580])
errorbar(wake_inds, PFC_1a_mean(wake_inds), PFC_1a_CI(1,wake_inds),PFC_1a_CI(2,wake_inds), ...
    'o', Color=[0.8500 0.3250 0.0980])
legend("sleep","wake");





%% Using the same matrices that hold all this data, plot histograms for CA1
% and PFC split into sleep and wake
CA1_1a_sleep = CA1_1a(:,sleep_inds); PFC_1a_sleep = PFC_1a(:,sleep_inds);
CA1_1a_srv = CA1_1a_sleep(~isnan(CA1_1a_sleep)); % sleep rate vector of all neurons.
PFC_1a_srv = PFC_1a_sleep(~isnan(PFC_1a_sleep));

% same for wake
CA1_1a_wake = CA1_1a(:,wake_inds); PFC_1a_wake = PFC_1a(:,wake_inds);
CA1_1a_wrv = CA1_1a_wake(~isnan(CA1_1a_wake)); 
PFC_1a_wrv = PFC_1a_wake(~isnan(PFC_1a_wake));

% Define bin edges for histograms, using the number of bins to be the sqrt
% of the smaller number of neurons between the wake and sleep sessions.
CA1_edges = linspace(0, max(CA1_1a,[],'all'), 3*sqrt(min([length(CA1_1a_srv),length(CA1_1a_wrv)])));
PFC_edges = linspace(0, max(PFC_1a,[],'all'), 3*sqrt(min([length(PFC_1a_srv),length(PFC_1a_wrv)])));

% KS test between the sleep and wake distributions
[~,pval_CA1] = kstest2(CA1_1a_srv,CA1_1a_wrv);
[~,pval_PFC] = kstest2(PFC_1a_srv,PFC_1a_wrv);

% Now plot histograms
figure;
sgtitle("Firing Rate Histograms")
subplot(2,2,1)
hold on;
title("CA1")
xlabel("Firing Rate (Hz)")
ylabel("Counts")
h_CA1s = histogram(CA1_1a_srv, CA1_edges, FaceColor=[0 0.01 1]);
h_CA1w = histogram(CA1_1a_wrv, CA1_edges, FaceColor=[0 0.5 0.8]);
txt = sprintf("pval = %d",pval_CA1);
text(max(CA1_edges)/4,max(h_CA1s.BinCounts)/2, txt)
legend("sleep","wake")

subplot(2,2,2)
hold on;
title("PFC")
xlabel("Firing Rate (Hz)")
ylabel("Counts")
h_PFCs = histogram(PFC_1a_srv, CA1_edges, FaceColor=[0.9 0.4 0]);
h_PFCw = histogram(PFC_1a_wrv, CA1_edges, FaceColor=[1 0.8 0]);
txt = sprintf("pval = %d",pval_PFC);
text(max(PFC_edges)/4,max(h_PFCs.BinCounts)/2, txt)
legend("sleep","wake")

% Plot 4-tile sublot of combinations of CA1 and PFC in sleep and wake
[~,pval_sleep] = kstest2(CA1_1a_srv,PFC_1a_srv);
[~,pval_wake] = kstest2(CA1_1a_wrv,PFC_1a_wrv);

subplot(2,2,3)
hold on;
title("Sleep")
xlabel("Firing Rate (Hz)")
ylabel("Counts")
h_CA1 = histogram(CA1_1a_srv, CA1_edges, FaceColor=[0 0.4470 0.7410]);
h_PFC = histogram(PFC_1a_srv, CA1_edges, FaceColor=[0.8500 0.3250 0.0980]);
txt = sprintf("pval = %d",pval_sleep);
text(max(CA1_edges)/4,max(h_CA1.BinCounts)/2, txt)
legend("CA1","PFC")

subplot(2,2,4)
hold on;
title("Wake")
xlabel("Firing Rate (Hz)")
ylabel("Counts")
h_CA1 = histogram(CA1_1a_wrv, CA1_edges, FaceColor=[0 0.4470 0.7410]);
h_PFC = histogram(PFC_1a_wrv, CA1_edges, FaceColor=[0.8500 0.3250 0.0980]);
txt = sprintf("pval = %d",pval_wake);
text(max(CA1_edges)/4,max(h_CA1.BinCounts)/2, txt)
legend("CA1","PFC")



%% Cumulative histograms

% Shantanu wants me to plot the histograms above, but in cumulative
% distributions.

% KS test between the CDFs. It turns out that the KS test of two samples is
% computed by turning them into CDFs and measuring the maximal difference
% between the two CDFs. So, I just need to hand kstest2 the data.
[~,pval_CA1] = kstest2(CA1_1a_srv,CA1_1a_wrv);
[~,pval_PFC] = kstest2(PFC_1a_srv,PFC_1a_wrv);
[~,pval_sleep] = kstest2(CA1_1a_srv,PFC_1a_srv);
[~,pval_wake] = kstest2(CA1_1a_wrv,PFC_1a_wrv);

figure;
sgtitle("Firing Rate CDFs")

subplot(2,2,1)
hold on;
CA1s_cdf = cdfplot(CA1_1a_srv);
CA1s_cdf.Color = [0 0.01 1];
CA1w_cdf = cdfplot(CA1_1a_wrv);
CA1w_cdf.Color=[0 0.5 0.8];
title("CA1")
ylabel("Fraction of Neurons")
xlabel("Firing Rate (Hz)")
legend("sleep","wake",Location='best')
txt = sprintf("pval = %.3d",pval_CA1);
text(CA1s_cdf.XData(end-1)*(2/4), 1/4, txt) % cdf line.XData has +- inf at
% the ends of the array, so to get the maximum plotted value index to
% end-1.

subplot(2,2,2)
hold on;
PFCs_cdf = cdfplot(PFC_1a_srv);
PFCs_cdf.Color = [0.9 0.4 0];
PFCw_cdf = cdfplot(PFC_1a_wrv);
PFCw_cdf.Color=[1 0.8 0];
title("PFC")
ylabel("Fraction of Neurons")
xlabel("Firing Rate (Hz)")
legend("sleep","wake",Location='best')
txt = sprintf("pval = %.3d",pval_PFC);
text(PFCs_cdf.XData(end-1)*(2/4), 1/4, txt);

subplot(2,2,3)
hold on;
CA1s_cdf = cdfplot(CA1_1a_srv);
CA1s_cdf.Color = [0 0.4470 0.7410];
PFCs_cdf = cdfplot(PFC_1a_srv);
PFCs_cdf.Color=[0.8500 0.3250 0.0980];
title("Sleep")
ylabel("Fraction of Neurons")
xlabel("Firing Rate (Hz)")
legend("CA1","PFC",Location='best')
txt = sprintf("pval = %.3d",pval_sleep);
text(CA1s_cdf.XData(end-1)*(2/4), 1/4, txt);

subplot(2,2,4)
hold on;
CA1w_cdf = cdfplot(CA1_1a_wrv);
CA1w_cdf.Color = [0 0.4470 0.7410];
PFCw_cdf = cdfplot(PFC_1a_wrv);
PFCw_cdf.Color=[0.8500 0.3250 0.0980];
title("Wake")
ylabel("Fraction of Neurons")
xlabel("Firing Rate (Hz)")
legend("CA1","PFC",Location='best')
txt = sprintf("pval = %.3d",pval_wake);
text(CA1w_cdf.XData(end-1)*(2/4), 1/4, txt);


%% Now I want to split the firing rates into the groups that Aanchal did:
% On track, wake during sleep sessions, NREM, and REM.




