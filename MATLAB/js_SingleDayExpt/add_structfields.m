function [] = add_structfields()
%ADD_STRUCTFIELDS   add fields to structs within a given data file across
%all rats and epochs. Overwrites the existing saved animal file and
%replaces it with the new version.
% Michael Satchell 08/30/2023

% Load data
data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
load_rats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','pos01'
filetypes = {'linfields01','cellinfo','spikes01'};

C_alldata = {}; % Cell array to hold data for all rats. If multiple filetypes 
% are loaded, each row holds a different file type, ordered in the same
% order as the elements in filetypes.

disp("Loading new animal data... ")
for r = 1:length(load_rats)
    fprintf("Loading animal: %s \n",load_rats{r})

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

linf_idx = find(contains(filetypes,'linfields'));
if isempty(linf_idx)
    error("linfields data must be loaded to run this analysis.")
end


% Extract spike data
C_allspikes = C_alldata(spikes_idx,:);

% Cell info data so I can determine brain region of each nrn. I could do a
% separate loop here and add an 'area' field to the neurons in C_allspikes.
% I could also just treat C_allinfo like C_allspikes and flatten the cells
% across tetrodes into one array inside the main for-loop, then grab the
% area from C_allinfo. For now I am going to choose the first option.
C_allinfo = C_alldata(cellinfo_idx,:);

% Linfield data
C_alllinf = C_alldata(linf_idx,:);



%% Add 'area', 'traj_x', and 'isPC' fields to spikes file.
% Add 'area' field to the spike data for later sorting. Also add place
% field data such as mean and std FR, sparsity, and spatial coverage for each trajectory.
% linfields must be loaded to add the placefield data.
for r = 1:size(C_allspikes,2)
    for e = 1:size(C_allspikes{1,r},2)
        for tet = 1:size(C_allspikes{1,r}{1,e},2)
            if ~isempty(C_allspikes{1,r}{1,e}{1,tet})
                for nrn = 1:size(C_allinfo{1,r}{1,e}{1,tet},2)
                    if ~isempty(C_allspikes{1,r}{1,e}{1,tet}{1,nrn}) 

                        % If the neuron exists in the spike file it should
                        % exist in the cellinfo file.
                        nrnCellinfo = C_allinfo{1,r}{1,e}{1,tet}{1,nrn};
                        

                        if isfield(nrnCellinfo,'area')
                            C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.area = nrnCellinfo.area; % Create field in spike struct.       
                        end
                        
                        % Linfields data is only available for behavioral
                        % epochs (even). Existence of neurons does not
                        % quite match that of C_allspikes.
                        if size(C_alllinf{1,r},2) >= e && ~isempty(C_alllinf{1,r}{1,e})
                            if ~isempty(C_alllinf{1,r}{1,e}{1,tet}) && size(C_alllinf{1,r}{1,e}{1,tet},2) >= nrn
                                if ~isempty(C_alllinf{1,r}{1,e}{1,tet}{1,nrn})
                                    trajsLinf = C_alllinf{1,r}{1,e}{1,tet}{1,nrn}; % Trajectories for a single neuron.
                                    for tr = 1:size(trajsLinf,2)
        
                                        trData = trajsLinf{tr};
                                        binRate = trData(:,7)./trData(:,6); % Firing rate normalized by occupancy
                                        FRmean = mean(binRate,1,'omitnan');
                                        FRstd = std(binRate,1,'omitnan');
        
                                        % calc sparsity
                                        totOcc = sum(trData(:,6),1); % Total time spent on track
                                        probOcc = trData(:,6)./totOcc; % Probability of occupying each spot on the trajectory
                                        sum1 = sum(probOcc.*binRate);
                                        sum2 = sum(probOcc.*(binRate.^2));
                                        sparsity = (sum1^2)/sum2;
        
                                        % calc spatial coverage
                                        numBins = size(trData,1); % number of spatial bins this trajectory is split into.
                                        peakRate = max(binRate);
                                        isAboveThr = binRate > 0.25*peakRate;
                                        coverage = sum(isAboveThr)/numBins; % Fraction of spatial bins with at least 25%
                                        % of the peak firing rate on that trajectory.

                                        
        
                                        % Assign to fields in struct
                                        C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_FRmeans(tr) = FRmean;
                                        C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_FRstds(tr) = FRstd;
                                        C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_sparsity(tr) = sparsity;
                                        C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_coverage(tr) = coverage;
                                        C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_FRpeaks(tr) = peakRate;
        
                                    end
                                    spStruct = C_allspikes{1,r}{1,e}{1,tet}{1,nrn};
                                    isPC = is_placecell(spStruct.traj_FRmeans, spStruct.traj_FRstds, spStruct.traj_FRpeaks,...
                                        spStruct.traj_sparsity, spStruct.traj_coverage);
                                    C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.isPlaceCell = isPC;

                                end
                            end
                        end

                    end
                end
            end
        end
    end

    % Save data from each rat to overwrite original file.
    spikes = C_allspikes(1,r);
    % save(string(fullfile(File_dir(f).folder, File_dir(f).name)), "spikes");


end