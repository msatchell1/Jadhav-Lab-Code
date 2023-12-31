% This script will recreate some results that Aanchal got on looking at
% firing rates of neurons in CA1 and PFC during awake and rest. It will
% also explore firing rate statistics in general.
%
% Michael Satchell 07/03/2023


%% Load data
clearvars;

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
load_rats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01'.
filetypes = {'cellinfo','spikes01','sleep01','waking01','sws01','rem01'};

C_alldata = {}; % Cell array to hold data for all rats. If multiple filetypes 
% are loaded, each row holds a different file type, ordered in the same
% order as the elements in filetypes.

disp("Loading new animal data... ")
for r = 1:length(load_rats)
    fprintf("Loaded animal: %s \n",load_rats{r})

    for ft = 1:length(filetypes)    

        short_name = load_rats{r};
        chop_idx = strfind(load_rats{r},'_') - 1;
        if ~isempty(chop_idx)
            short_name = load_rats{r}(1:chop_idx); % Gets the first characters of the rat's name before an '_'.
            % So far this is only needed for ER1_NEW to remove the '_NEW'.
        end

        % Does not load from EEG folder
        File_dir = dir(data_dir+"/"+load_rats(r)+'_direct'+"/"+short_name+filetypes{ft}+"*");
    
        if isempty(File_dir)
            error("%s file does not exist for animal: %s \n",filetypes{ft},load_rats{r})
        elseif length(File_dir) > 1
            error("More than one file detected when searching for: %s, in animal: %s \n" + ...
                "Change names in filetypes to be more specific.", filetypes{ft},load_rats{r});
        else
            file = struct2cell(load(string(fullfile(File_dir.folder, File_dir.name)))); % load data
            file = file{:};
            C_alldata{ft,r} = file{1,1};
            fprintf("       Loaded file: %s  \n",  File_dir.name)
        end
    end
end

C_alldata = clip_17_epochs(C_alldata); % removes extra epoch data.



%% Plot mean and median of mean firing rates for one rat
% NOTE THE CODE BELOW ONLY WORKS WHEN 'cellinfo' IS THE FIRST ITEM IN
% filetypes.

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
for a = 1:size(C_alldata,2)
    epochcounts(a) = length(C_alldata{1,a});
end

% First I need to assemble the firing rate information from all the rats
% into one matrix for CA1 and one for PFC.
CA1_allrates = {}; PFC_allrates = {};

% Loop all animals
for a = 1:size(C_alldata,2)


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

CI_val = 75; % Confidence interval to be used in %.

for e = 1:size(CA1_1a_mean,2)
    % CI from a normal? distribution with N-1 data points
    temp_CI_CA1 = tinv([(1-CI_val/100)/2, 1-(1-CI_val/100)/2], CA1_numnrns(1,e)-1); 
    temp_CI_PFC = tinv([(1-CI_val/100)/2, 1-(1-CI_val/100)/2], PFC_numnrns(1,e)-1);

    % multiplying this CI by the SEM to get the error size
    CA1_1a_CI(:,e) = bsxfun(@times, CA1_1a_SEM(e), temp_CI_CA1(:));
    PFC_1a_CI(:,e) = bsxfun(@times, PFC_1a_SEM(e), temp_CI_PFC(:));
end

% Indices for rest and wake
rest_inds = 1:2:size(CA1_1a_mean,2);
wake_inds = 2:2:size(CA1_1a_mean,2);

figure;
sgtitle(sprintf("Mean Firing Rates Across All Neurons From %d Animals",length(load_rats)));

subplot(1,2,1)
hold on;
title("CA1")
ylabel("Mean Firing Rate (Hz)")
xlabel("Epoch")
errorbar(rest_inds, CA1_1a_mean(rest_inds), CA1_1a_CI(1,rest_inds),CA1_1a_CI(2,rest_inds), ...
    'o',Color=[0 0.4470 0.7410], MarkerFaceColor='blue')
errorbar(wake_inds, CA1_1a_mean(wake_inds), CA1_1a_CI(1,wake_inds), CA1_1a_CI(2,wake_inds), ...
    'o', Color=[0 0.4470 0.7410])
legend("rest","wake");

subplot(1,2,2)
hold on;
title("PFC")
ylabel("Mean Firing Rate (Hz)")
xlabel("Epoch")
errorbar(rest_inds, PFC_1a_mean(rest_inds), PFC_1a_CI(1,rest_inds),PFC_1a_CI(2,rest_inds), ...
    'o',Color=[0.8500 0.3250 0.0980], MarkerFaceColor=[0.8500 0.3250 0.0580])
errorbar(wake_inds, PFC_1a_mean(wake_inds), PFC_1a_CI(1,wake_inds),PFC_1a_CI(2,wake_inds), ...
    'o', Color=[0.8500 0.3250 0.0980])
legend("rest","wake");





% Using the same matrices that hold all this data, plot histograms for CA1
% and PFC split into rest and wake. 
CA1_1a_rest = CA1_1a(:,rest_inds); PFC_1a_rest = PFC_1a(:,rest_inds);
CA1_1a_rrv = CA1_1a_rest(~isnan(CA1_1a_rest)); % rest rate vector of all neurons.
PFC_1a_rrv = PFC_1a_rest(~isnan(PFC_1a_rest));

% same for wake
CA1_1a_wake = CA1_1a(:,wake_inds); PFC_1a_wake = PFC_1a(:,wake_inds);
CA1_1a_wrv = CA1_1a_wake(~isnan(CA1_1a_wake)); 
PFC_1a_wrv = PFC_1a_wake(~isnan(PFC_1a_wake));

% Define bin edges for histograms, using the number of bins to be the sqrt
% of the smaller number of neurons between the wake and rest sessions.
CA1_edges = linspace(0, max(CA1_1a,[],'all'), 3*sqrt(min([length(CA1_1a_rrv),length(CA1_1a_wrv)])));
PFC_edges = linspace(0, max(PFC_1a,[],'all'), 3*sqrt(min([length(PFC_1a_rrv),length(PFC_1a_wrv)])));

% KS test between the rest and wake distributions
[~,pval_CA1] = kstest2(CA1_1a_rrv,CA1_1a_wrv);
[~,pval_PFC] = kstest2(PFC_1a_rrv,PFC_1a_wrv);

% Now plot histograms
figure;
sgtitle("Firing Rate Histograms")
subplot(2,2,1)
hold on;
title("CA1")
xlabel("Firing Rate (Hz)")
ylabel("Counts")
h_CA1s = histogram(CA1_1a_rrv, CA1_edges, FaceColor=[0 0.01 1]);
h_CA1w = histogram(CA1_1a_wrv, CA1_edges, FaceColor=[0 0.5 0.8]);
txt = sprintf("pval = %d",pval_CA1);
text(max(CA1_edges)/4,max(h_CA1s.BinCounts)/2, txt)
legend("rest","wake")

subplot(2,2,2)
hold on;
title("PFC")
xlabel("Firing Rate (Hz)")
ylabel("Counts")
h_PFCs = histogram(PFC_1a_rrv, CA1_edges, FaceColor=[0.9 0.4 0]);
h_PFCw = histogram(PFC_1a_wrv, CA1_edges, FaceColor=[1 0.8 0]);
txt = sprintf("pval = %d",pval_PFC);
text(max(PFC_edges)/4,max(h_PFCs.BinCounts)/2, txt)
legend("rest","wake")

% Plot 4-tile sublot of combinations of CA1 and PFC in rest and wake
[~,pval_rest] = kstest2(CA1_1a_rrv,PFC_1a_rrv);
[~,pval_wake] = kstest2(CA1_1a_wrv,PFC_1a_wrv);

subplot(2,2,3)
hold on;
title("Rest")
xlabel("Firing Rate (Hz)")
ylabel("Counts")
h_CA1 = histogram(CA1_1a_rrv, CA1_edges, FaceColor=[0 0.4470 0.7410]);
h_PFC = histogram(PFC_1a_rrv, CA1_edges, FaceColor=[0.8500 0.3250 0.0980]);
txt = sprintf("pval = %d",pval_rest);
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



% Cumulative histograms

% Shantanu wants me to plot the histograms above, but in cumulative
% distributions.

% KS test between the CDFs. It turns out that the KS test of two samples is
% computed by turning them into CDFs and measuring the maximal difference
% between the two CDFs. So, I just need to hand kstest2 the data.
[~,pval_CA1] = kstest2(CA1_1a_rrv,CA1_1a_wrv);
[~,pval_PFC] = kstest2(PFC_1a_rrv,PFC_1a_wrv);
[~,pval_rest] = kstest2(CA1_1a_rrv,PFC_1a_rrv);
[~,pval_wake] = kstest2(CA1_1a_wrv,PFC_1a_wrv);

figure;
sgtitle("Firing Rate CDFs")

subplot(2,2,1)
hold on;
CA1s_cdf = cdfplot(CA1_1a_rrv);
CA1s_cdf.Color = [0 0.01 1];
CA1w_cdf = cdfplot(CA1_1a_wrv);
CA1w_cdf.Color=[0 0.5 0.8];
title("CA1")
ylabel("Fraction of Neurons")
xlabel("Firing Rate (Hz)")
legend("rest","wake",Location='best')
txt = sprintf("pval = %.3d",pval_CA1);
text(CA1s_cdf.XData(end-1)*(2/4), 1/4, txt) % cdf line.XData has +- inf at
% the ends of the array, so to get the maximum plotted value index to
% end-1.
set(gca,'XScale','log');

subplot(2,2,2)
hold on;
PFCs_cdf = cdfplot(PFC_1a_rrv);
PFCs_cdf.Color = [0.9 0.4 0];
PFCw_cdf = cdfplot(PFC_1a_wrv);
PFCw_cdf.Color=[1 0.8 0];
title("PFC")
ylabel("Fraction of Neurons")
xlabel("Firing Rate (Hz)")
legend("rest","wake",Location='best')
txt = sprintf("pval = %.3d",pval_PFC);
text(PFCs_cdf.XData(end-1)*(2/4), 1/4, txt);
set(gca,'XScale','log');

subplot(2,2,3)
hold on;
CA1s_cdf = cdfplot(CA1_1a_rrv);
CA1s_cdf.Color = [0 0.4470 0.7410];
PFCs_cdf = cdfplot(PFC_1a_rrv);
PFCs_cdf.Color=[0.8500 0.3250 0.0980];
title("Rest")
ylabel("Fraction of Neurons")
xlabel("Firing Rate (Hz)")
legend("CA1","PFC",Location='best')
txt = sprintf("pval = %.3d",pval_rest);
text(CA1s_cdf.XData(end-1)*(2/4), 1/4, txt);
set(gca,'XScale','log');

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
set(gca,'XScale','log');


%% Split the firing rates into the following states during rest epochs:
% Behavior(on track), wake during sleep sessions, sleep, NREM, and REM.

% Because cellinfo just holds the mean rate for a cell during an entire
% epoch, I can't use this to evaluate rest states. Instead I need to use
% the spiking data for each cell, split it by rest state, and take the
% mean during all instances of that state in an epoch. What I want is to
% come up with a matrix like that used above, where each column is an epoch
% and cell spiking data is held in the rows. It would be nice to set it up
% so that each row represents a cell that is consistent across epochs. 

% I did a check: the number of tetrodes used to detect spikes varies
% between animals, even between the 32-tetrode ones. This is likely because
% some tetrodes stopped working, or were used for other purposes. The
% actual number of tetrodes in data files like spikes01 varies between about 29
% and 32 for the 32-tetrode animals. 
% Also, thankfully the data processing pipeline keeps a consistent number
% of cells listed on each tetrode across epochs. Even though some neurons
% are lost or gained over the course of the epochs, those neurons just have
% an empty array as a placeholder. This means that cell identity is
% maintained as the column under each tetrode across epochs. This is useful
% because now I can know the size of the matrices I need to make with
% number of cells as the first dimension, and this will be consistent
% across epochs, so I can use epochs as the second dimension.

spikes_idx = find(contains(filetypes,'spikes01')); % Gets row index of spike data in C_alldata.
% Note that indexing like this is only supposed to give one index as a
% result.
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end

cellinfo_idx = find(contains(filetypes,'cellinfo'));
if isempty(cellinfo_idx)
    error("cellinfo data must be loaded to run this analysis.")
end

% List of all states to sort the spike data into. Note that adding
% rippletime01 adds quite a bit of calculation time. 
% WARNING: This does not change the order of how these states are stored in
% C_allstates or later FR_allStates. Thus, stateNames MUST MATCH ORDER
% THESE FILES ARE LISTED IN filetypes.
stateNames = {'sleep','waking','sws','rem'};

states_idx = find(contains(filetypes, stateNames));

% for s = 1:length(stateNames)
%     states_idx(s) = find(contains(filetypes, stateNames{s}));
% end
if isempty(states_idx)
    error("at least one of the following rest state data must be loaded: \n +" + ...
        "   sleep01, waking01, sws01, rem01, or rippletime01.")
end

brainAreas = {'CA1','PFC'}; % Brain areas to split data into. Should be only CA1 and PFC

% Extract spike data
C_allspikes = C_alldata(spikes_idx,:);

% Cell info data so I can determine brain region of each nrn. I could do a
% separate loop here and add an 'area' field to the neurons in C_allspikes.
% I could also just treat C_allinfo like C_allspikes and flatten the cells
% across tetrodes into one array inside the main for-loop, then grab the
% area from C_allinfo. For now I am going to choose the first option.
C_allinfo = C_alldata(cellinfo_idx,:);

% Add 'area' field to the spike data for later sorting.
for r = 1:length(C_allspikes)
    for e = 1:length(C_allspikes{1,r})
        for tet = 1:length(C_allspikes{1,r}{1,e})
            if ~isempty(C_allspikes{1,r}{1,e}{1,tet})
                for nrn = 1:length(C_allinfo{1,r}{1,e}{1,tet})
                    if ~isempty(C_allspikes{1,r}{1,e}{1,tet}{1,nrn})

                        % If the cell exists in the spike file it should
                        % exist in the cellinfo file.
                        nrnCellinfo = C_allinfo{1,r}{1,e}{1,tet}{1,nrn};
                        if isfield(nrnCellinfo,'area')
                            C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.area = nrnCellinfo.area; % Create field in spike struct.
                           
                        end
                    end
                end
            end
        end
    end
end

% Store all state data in one array. The order of the data
% stored in rows is that of stateNames.
C_allstates = C_alldata(states_idx,:);


% Cell array to hold matrices for each brain region (in brainArea order).
% Matrices contain mean spike rates by epoch, not discriminating by state. 
% Rows are neuron identity, columns are epochs. This is for all rats.
FR_byEpoch = cell(1,length(brainAreas));

for r = 1:length(C_allspikes) % Nested loops to get to each neuron.
    
    for a = 1:length(brainAreas) % I could eliminate this loop by using strfind(brainAreas) to index...

        % Creates a matrix to hold spike rate data for one rat. Dims (num nrns)
        % x (num epochs). The first epoch is used to find the number of nrns
        % because this does not change across epochs as stated above. Third
        % dimension is area in order brainAreas.
        ratRates = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));

        for e = 1:length(C_allspikes{1,r})
    
            nrnsAlltets = [C_allspikes{1,r}{1,e}{:}]; % Combining nrn data from all tets.
    
            for nrn = 1:length(nrnsAlltets)
    
                % isfield returns false if the struct does not exist
                if isfield(nrnsAlltets{nrn},'data')
                    
                    if strcmp(nrnsAlltets{nrn}.area, brainAreas{a})
                        % I checked and confirmed that calculating the mean firing rate 
                        % this way matches the .meanrate field value for each nrn (only off by 0.0001 for some nrns).
                        ratRates(nrn,e) = size(nrnsAlltets{nrn}.data,1)/(nrnsAlltets{nrn}.timerange(2)-nrnsAlltets{nrn}.timerange(1));
                    end
                end
            end
        end
        FR_byEpoch{a} = [FR_byEpoch{a}; ratRates]; 
    end
end

% Separating into rest and behavior epoch data.
behEpochs = 2:2:size(FR_byEpoch{1},2);
restEpochs = 1:2:size(FR_byEpoch{1},2);
FRbeh = {}; FRrest = {};
for a = 1:size(FR_byEpoch,2)
    FRbeh{a} = FR_byEpoch{a}(:,behEpochs);
    FRrest{a} = FR_byEpoch{a}(:,restEpochs);
end


% Gets epochs that have data for all states. For most states this should be all odd
% epochs for the data I am sorting here, because all the states
% happen during the rest epochs except for rippletimes. I will detect non-empty
% epochs for generality. Note this assumes all rats have the same number of
% epochs as rat 1, which should be true when clip_17_epochs.m has been
% used on C_allstates.
nonemptyEpochs = zeros(size(C_allstates,1),length(C_allstates{1,1}));
for s = 1:size(C_allstates,1)
    nonemptyEpochs(s,:) = ~cellfun(@isempty, C_allstates{s,1});
end


% To hold firing rates for all states held in C_allstates
FR_allStates = cell(size(C_allstates,1),length(brainAreas));
% Time duration of states by epoch. (num rats) x (num states) x (num epochs).
Dur_allStates = zeros(size(C_allstates,2),size(C_allstates,1),size(nonemptyEpochs,2)); 
% Main loop to calculate mean firing rate for all states in all brain areas.
for a = 1:length(brainAreas)

    for s = 1:size(C_allstates,1)
    
        % Gets epochs that have data for this state. This should be identical across rats.
        nonempty = zeros(1,length(C_allstates{1,1}));
    
        for r = 1:size(C_allstates,2) 
            % Indicates which epochs contain data
            nonempty(1,:) = ~cellfun(@isempty, C_allstates{s,r});
            % To hold mean firing rate data across all epochs for a single rat.
            ratMeans = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));
    
            for e = 1:length(C_allstates{s,r})
                
                if nonempty(1,e) == 1
                % Now the real deal. For every non-empty epoch, each state type (rem, rippletime,
                % sleep, sws, and waking) will loop through the times
                % that that state occurs, counting the number of spikes
                % that occur during the duration of that occurance. In the end
                % these spike counts will be summed up and divided by the
                % total duration of that state to get the mean rate.
                
                epochData = C_allstates{s,r}{1,e}; % Data on a single rat, state, and epoch.
                nrnsAlltets = [C_allspikes{1,r}{1,e}{:}]; % Combining nrn spike data from all tets.
                STs_inState = cell(1,length(nrnsAlltets)); % Cell of spike times that occur within the given state.
                % Columns are neurons.
    
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
                
                stateDur = sum(epochData.endtime - epochData.starttime);  % Total time spent in state.
                Dur_allStates(r,s,e) = stateDur;
                FRmeans = NaN(size(STs_inState)); % To hold mean firing rate of nrns.
                for nrn = 1:size(STs_inState,2)
                    if ~isempty(STs_inState{1,nrn}) % Important to leave nrns with no spikes as NaNs.
                        meanFR = length(STs_inState{1,nrn})/stateDur; % Calculate mean firing rate
                        FRmeans(1,nrn) =  meanFR;
                    end
                end
                 
                ratMeans(:,e) = FRmeans; % Stores all FR means for that epoch.
               
                end
            end
            FR_allStates{s,a} = [FR_allStates{s,a}; ratMeans];   
            
        end
    end
end





%% Plot the rates as distributions
stateHistF = figure;
sgtitle("Mean FR Distributions During Rest Epochs")
faceColors = [];

for s = 1:size(FR_allStates,1)
   
    edges = linspace(0, max(FR_allStates{s},[],'all'), 10*sqrt(sum(~isnan(FR_allStates{s}),'all')));
    subplot(ceil(size(FR_allStates,1)/2),2,s)
    hold on;
    hCA1 = histogram(FR_allStates{s,1}(~isnan(FR_allStates{s,1})), edges);
    hPFC = histogram(FR_allStates{s,2}(~isnan(FR_allStates{s,2})), edges);
    title(stateNames{s})
    ylabel("Counts")
    xlabel("Firing Rate (Hz)")
    legend(brainAreas)
end
linkaxes(findobj(stateHistF,'Type','axes'), 'x');



% Plot the CDF for each brain region together
stateCDFF = figure;
sgtitle("Cumulative Histogram Curves Separated by Rest State (Rest Epochs Only)")
colors = {[0 0.4470 0.7410],[1 0.9 0],[0.5 0.1 1],[1 0 0]};

for a = 1:size(FR_allStates,2)

    subplot(1,size(FR_allStates,2),a)
    hold on;
    % Plot all rest states
    for s = 1:size(FR_allStates,1)
        cdf = cdfplot(FR_allStates{s,a}(~isnan(FR_allStates{s,a})));
        cdf.Color = colors{s}; cdf.LineWidth=2;
    end
    % Plot behavior
    cdf = cdfplot(FRbeh{1,a}(~isnan(FRbeh{1,a})));
    cdf.Color = [0.3 0.9 0]; cdf.LineWidth=2;

    set(gca,'XScale','log');
    title(brainAreas{a})
    xlabel("Firing Rate (Hz)")
    ylabel("Fraction of Neurons")
    legend([stateNames,"beh"],Location='best')

end


% Histograms of all data from behavior or rest epochs.
figure;
sgtitle("Mean FR Distributions by Epoch Type")

axbeh = subplot(1,2,1);
hold on;
edges = linspace(0, max(FRbeh{1},[],'all'), 10*sqrt(sum(~isnan(FRbeh{1}),'all')));
histogram(FRbeh{1}(~isnan(FRbeh{1})), edges);
histogram(FRbeh{2}(~isnan(FRbeh{2})), edges);
title("Behavior (even epochs)")
ylabel("Counts")
xlabel("Firing Rate (Hz)")
legend(brainAreas)

axrest = subplot(1,2,2);
hold on;
edges = linspace(0, max(FRrest{1},[],'all'), 10*sqrt(sum(~isnan(FRrest{1}),'all')));
histogram(FRrest{1}(~isnan(FRrest{1})), edges);
histogram(FRrest{2}(~isnan(FRrest{2})), edges);
title("Rest (odd epochs)")
ylabel("Counts")
xlabel("Firing Rate (Hz)")
legend(brainAreas)

linkaxes([axbeh,axrest],'x');





% Plot rates by epoch and state
stateEpochF = figure;
sgtitle("Mean FR Across Epochs")
colors = {[0 0.4470 0.7410],[1 0.9 0],[0.5 0.1 1],[1 0 0]}; % Color of data, will be plotted with stateNames order.
for a = 1:size(FR_allStates,2)
   
    subplot(1,size(FR_allStates,2),a)
    hold on;

    for s = 1:size(FR_allStates,1)
        % Calculate mean firing rate, standard deviation, and SEM values for
        % each rest epoch.
        FRmeans = mean(FR_allStates{s,a},1,'omitnan');
        FRstds = std(FR_allStates{s,a},0,1,'omitnan');
        numnrns = sum(~isnan(FR_allStates{s,a}),1);
        FRsems = FRstds./sqrt(numnrns);
        
        FR_CIs = zeros(2,size(FRmeans,2));
        CI_pct = 75; % Confidence interval value in percent.
        for e = 1:size(FRmeans,2)
            temp_CI = tinv([(1-CI_pct/100)/2, 1-(1-CI_pct/100)/2], numnrns(1,e)-1); 
            FR_CIs(:,e) = bsxfun(@times, FRsems(e), temp_CI(:));
        end

        errorbar(restEpochs,FRmeans(restEpochs),FR_CIs(1,restEpochs),FR_CIs(2,restEpochs),...
            ".", Color=colors{s})
        plot(restEpochs,FRmeans(restEpochs), Color=colors{s}, LineWidth=2, HandleVisibility='off')
    end
    
    % Doing the same for wake epochs.
    FRmeans = mean(FRbeh{1,a},1,'omitnan');
    FRstds = std(FRbeh{1,a},0,1,'omitnan');
    numnrns = sum(~isnan(FRbeh{1,a}),1);
    FRsems = FRstds./sqrt(numnrns);
    FR_CIs = zeros(2,size(FRmeans,2));
    for e = 1:size(FRmeans,2)
            temp_CI = tinv([(1-CI_pct/100)/2, 1-(1-CI_pct/100)/2], numnrns(1,e)-1); 
            FR_CIs(:,e) = bsxfun(@times, FRsems(e), temp_CI(:));
    end

    errorbar(behEpochs,FRmeans,FR_CIs(1,:),FR_CIs(2,:),...
            ".", Color=[0.3 0.9 0])
    plot(behEpochs,FRmeans, Color=[0.3 0.9 0], LineWidth=2, HandleVisibility='off')

    title(brainAreas{a})
    ylabel("Mean Firing Rate (Hz)")
    xlabel("Epoch")
    L = legend([stateNames,"beh"]);

end


% Plot time duration of states across epochs
figure;
colors = {[0 0.4470 0.7410],[1 0.9 0],[0.5 0.1 1],[1 0 0]}; % Color of data, will be plotted with stateNames order.
hold on;

DurMean = squeeze(mean(Dur_allStates,1,'omitnan')); % Take the mean across all rats. Remove dims of length 1
DurStd = squeeze(std(Dur_allStates,0,1,'omitnan'));

for s = 1:size(DurMean,1)
    errorbar(restEpochs,DurMean(s,restEpochs)./60, DurStd(s,restEpochs)./60, '.',Color=colors{s})
    plot(restEpochs,DurMean(s,restEpochs)./60,Color=colors{s}, LineWidth=2, HandleVisibility='off')
end
title("Duration of States Across Animals (Rest Epochs Only)")
ylabel("Time in State (min)")
xlabel("Epoch")
legend(stateNames)
