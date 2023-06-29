% This script is just for exploring and getting familiar with the data from
% Justin's single day W track experiment.
%
% Michael Satchell 6/26/23


%% Notes on data format and content

% The good animals to do data analysis on are: BG1, ER1_NEW, JS
% 12,13,14,15,17,21,34, KC10, KL8, and ZT2.

% JS20 does not have many cells but behaved well. KC1 did not learn the
% task. Both of these animals do not have _direct files for them, meaning
% that the preprocessing of the data has not been done. 

% The format for most of the data files is as follows:
% Day / Epoch / Tetrode / Cell 

% There are 17 epochs, with 8 W track behavior sessions and 9 sleep
% sessions. Every odd epoch is a sleep epoch. Most animals have 32 tetrodes,
% although I think at least one has 64 tetrodes.


% Each animal has an EEG folder which inside contains a ton of files. These
% files are under the data type categories:
% delta, eeg, eegref, gamma, ripple, theta
% Under each of these catagories, there are num_epochs x num_tetrodes
% files. So for the 32 tetrode animals there are 544 files under each of
% the data type categories.
% The file names are formatted as such:
% (ANIMAL NAME)(data type)(day)-(epoch)-(tetrode).mat
% Since all the experiments were conducted in one day for each animal,
% (day) always reads 01 (or at least it should I think).
%
% All files under the EEG folder follow the following data format in the
% loaded .mat file:
% day / epoch / tetrode
% For each file, all cells are empty except that which corresponds to the 
% epoch/tetrode indicated in the filename. The struct inside the .mat file
% contains all the information for the file, with the field 'data'
% containing all the LFP data.
%
% EEG (raw LFP) is data recorded from the electrode with reference to a 
% ground screw. EEGref is the same data, but is recorded with respect 
% to a local reference: an electrode in a quiet area (like corpus collosum).
% In other words, both EEG and EEGref hold the data recorded from the 
% tetrode of interest, however, this data is recorded in reference to either
% the ground screw to eliminate environmental noise (EEG) or in reference 
% to the reference electrode in a quiet area of the brain to eliminate 
% similarities between the two brain areas (EEGref). This LFP trace with reference to an
% electrode in a quiet area of the brain is needed because ground screw 
% referencing does not remove some high frequency noise that is, say, only 
% detected on the tetrode in CA1 but not on the ground screw. The EEGref 
% trace will more often catch the high frequency noise, and so is better to 
% use when looking at high frequency signals such as ripples. For all low 
% frequency things I can just use EEG. EEGref can be used for more universal
% measures of oscillations like theta, because the oscillations measured 
% there are considered more global to the whole hippocampus, especially 
% when considering phase. For example, the theta phase of many tetrodes 
% might be misaligned, so how do we know what theta phase to use when measuring
% phase precession of a population? We can measure with respect to the theta
% on the reference electrode to get a more universal template of the oscillation.  
% 
% The delta, gamma, ripple,
% and theta files are filtered versions of the eeg files. Justin said that
% only the ripple data files need to have the reference LFP trace
% subtracted from it because noise is usually only high at that frequency.
%
% Ripples are calculated from tetrodes with cells, and are usually coherent
% across all tetrodes so it doesn't really matter which tetrode you use to
% measure the ripples. Ripples are also only detected for low velocity
% periods. I am assuming that the EEGref trace is used when finding ripples
% on a tetrode.


% Outside the EEG folder are a bunch of other files which contain more
% processed data, cell data, spike data, etc:
% cellinfo - Information about cells that were clustered (clustered meaning
% identified as cells using mountainsort)
% DIO - behavior triggers
% linfields - Main file for place fields. Includes place field metrics.
% task - Contains linearcoord which I can use to reference between the
% linearized coordinates and the 2D spatial map.
% sws - Contains start and end time of slow wave sleep. SWS is detected as
% follows: if the animal doesn't move for 1 min in the sleep cage, it is
% assumed to be asleep and entered into SWS. If the theta/delta ratio goes
% above a certain threshold then REM sleep is said to be occuring.
% rem - start and end times of REM sleep.
% spikes - Contains spike information for each cell. The column each neuron
% sits in under a tetrode in spikes matches that in cellinfo. And the cell
% identity is conserved across epochs.
% rippletime - contains start and end times for fully processed ripples. These are ripples
% that have been collapsed across the multiple ripple-detecting tetrodes in
% the CA1 pyramidal layer. This also excludes some ripple events based
% on rat velocity. Note the rippletime file is missing for some animals.


% More detailed information on the format of specific file types is as
% follows:
%
%
% linfields:
% day / epoch / tetrode / cell / trajectory
% There are four columns here representing the four W-track trajectories.
% Col 1: out right, col 2: in right, col 3: out left, col 4: in left.
% Under each trajectory is the place field data for that cell on that
% trajectory. Column 1 is the distance in cm along that linearized
% trajectory. Column 5 contains the occupancy normalized smooth firing field.
%
% LFP files:
% day / epoch / tetrode 
% 'samprate' is the rate at which the data in this file was recorded at.
% For the filter data files, under 'data' there will be 3 columns. Column 1
% is the amplitude of the filtered signal (in mV) multiplied by the scaling
% factor 'voltage_scaling'. To get back to mV, simply divide the amplitude
% by this factor. Column 2 holds instantaneous phase (in radians) multipled
% by 10,000. Column 3 holds the filter envelope magnitude.
%
% spikes:
% day / epoch / tetrode / cell
% 'data' holds the spiking data for that cell. Each row represents one
% spike for that cell. Column 1 is the time of that spike, column 2 the
% 2D x-position at the time of the spike, column 3 the 2D y-position, 
% column 4 is the animal's head direction in radians, column 5 is unused,
% column 6 is the largest amplitude of the spike measured by any of the
% four electrodes on the tetrode, column 7 is the linearized position
% index.
%
% cellinfo:
% day / epoch / tetrode / cell
% 'meanrate' is mean firing rate of that cell during that epoch.
% 'numspikes' is number of spikes during that epoch. 'csi' is the complex
% spike index, which shows how much the amplitude of spikes attenuates during 
% a burst. 'propbursts' are the proportion of spikes that are in bursts, with 
% bursts defined as spikes with < 6ms ISIs.
% Some cells during some epochs have all but the 'tag2' field missing, with this
% field containing 'Undef'...
% Justin says these are cells that have clusters in other epochs (by
% clusters he means being identifiable as a cell in mountainsort or
% matclust). In fact, he says that in CA1 there is a large shift in the
% population of active cells between wake and sleep, such that many cells
% are only active during waking epochs and others only during sleep. 




%% Load data

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
load_rats = {'BG1_direct'}; % Which rats to load data from
filetype = {'linfields'};
%filetype = {'ripple01-03-30', 'delta01-03-30', 'gamma01-03-30', 'theta01-03-30', 'eeg01-03-30', 'eegref01-03-30'}; % File type to load for each rat. All files will 
% be loaded that contain filetype in the file name. Use 'eeg01' when you
% want all the LFP files, otherwise eegref will also be included.



for i=1:length(load_rats)
    
    for k = 1:length(filetype)

        File_dir = dir(data_dir+"/"+load_rats(i)+"/**/*"+filetype(k)+"*");
        for j = 1:length(File_dir)
            load(string(fullfile(File_dir(j).folder, File_dir(j).name)));
        end

    end

end





%% Plotting for linfields

e2 = linfields{1,1}{1,2}; % Data from second epoch
e2tet2 = e2{1,2}; % Data from second tetrode
e2tet2cell1 = e2tet2{1,1}; % Data from first cell

traj1 = e2tet2cell1{1,1};

figure;
hold on;
title("Smooth, Linearized, Occupancy Normalized Place Field")
ylabel("Firing Rate (Hz)??")
xlabel("Linearized Position (cm)")
plot(traj1(:,1), traj1(:,5))




%% Plotting for LFP data

% Assuming LFP data from 01-03-30 was loaded:
S_eeg = eeg{1,1}{1,3}{1,30};
time_range = S_eeg.starttime:(1/S_eeg.samprate):S_eeg.endtime;
S_eegref = eegref{1,1}{1,3}{1,30};
S_delta = delta{1,1}{1,3}{1,30};
S_theta = theta{1,1}{1,3}{1,30};
S_ripple = ripple{1,1}{1,3}{1,30};
S_gamma = gamma{1,1}{1,3}{1,30};


% figure;
% hold on;
% title("Raw LFP Trace")
% ylabel("Voltage (mV)??")
% xlabel("Time (s)")
% plot(time_range,S_eeg.data)


figure;
hold on;
title("Some LFP Traces Together")
xlabel("Time (s)")
plot(time_range, S_eeg.data, DisplayName='raw LFP')
plot(time_range, S_eegref.data, DisplayName='ref LFP')
% plot(time_range, S_delta.data(:,1),DisplayName='delta')
legend();

% figure;
% hold on;
% title("Filtered LFP Traces")
% xlabel("Time (ms)")
% plot(time_range, S_delta.data(:,1),DisplayName='delta')
% plot(time_range, S_gamma.data(:,1),DisplayName='gamma')
% plot(time_range, S_ripple.data(:,1),DisplayName='ripple')
% plot(time_range, S_theta.data(:,1),DisplayName='theta')
% legend()

%% Subplot showing individual filters
figure;
subplot(2,2,1)
hold on;
title("Delta Filter (1-4 Hz)")
xlabel("Time (ms)")
plot(time_range, S_eeg.data, DisplayName='raw LFP')
plot(time_range, S_delta.data(:,1),DisplayName='delta')
legend()

subplot(2,2,2)
hold on;
title("Theta Filter (6-12 Hz)")
xlabel("Time (ms)")
plot(time_range, S_eeg.data, DisplayName='raw LFP')
plot(time_range, S_theta.data(:,1),DisplayName='theta')
legend()

subplot(2,2,3)
hold on;
title("Gamma Filter (40-100 Hz)")
xlabel("Time (ms)")
plot(time_range, S_eeg.data, DisplayName='raw LFP')
plot(time_range, S_gamma.data(:,1),DisplayName='gamma')
legend()

subplot(2,2,4)
hold on;
title("Ripple Filter (150-250 Hz)")
xlabel("Time (ms)")
plot(time_range, S_eeg.data, DisplayName='raw LFP')
plot(time_range, S_ripple.data(:,1),DisplayName='ripple')
legend()



%% Cell Spiking Data

load_rats = {'BG1_direct'}; % Which rats to load data from
filetype = {'spikes', 'cellinfo'};

for i=1:length(load_rats)

    for k = 1:length(filetype)

        File_dir = dir(data_dir+"/"+load_rats(i)+"/**/*"+filetype(k)+"*");
        for j = 1:length(File_dir)
            load(string(fullfile(File_dir(j).folder, File_dir(j).name)));
        end

    end

end


% % Grab data from epoch 3 (second sleep epoch), second tetrode, second cell
% S_spikes = spikes{1,1}{1,3}{1,2}{1,2};
% spikedata = S_spikes.data;
% 
% S_cell = cellinfo{1,1}{1,3}{1,2}{1,2};

% figure;
% hold on;
% scatter(spikedata(:,1),zeros(size(spikedata(:,1))),'red', ".") %, spikedata(:,6))



% Now say I want to make a raster plot of all cells in an epoch on an animal

e3spikes = spikes{1,1}{1,3};
e3spikes_unpacked = [e3spikes{:}];
S_e3spikes = [e3spikes_unpacked{:}];

e3info = cellinfo{1,1}{1,3};
e3info_unpacked_all = [e3info{:}];

infodataexist = [];
% Some of the structs don't have complete fields with data. To remove
% these, I will only retain structs that have the 'meanrate' field.
for i = 1:length(e3info_unpacked_all)
    infodataexist(i) = isfield(e3info_unpacked_all{i}, 'meanrate');
end

infodataexist = logical(infodataexist); % For some reason the values in the array are double 
% and not logical, even though isfield returns logical...
e3info_unpacked = e3info_unpacked_all(infodataexist);
% Removing the structs with no data makes the number of cells we have
% meaningful info on equal to the number of cellswith spike data. So
% cellinfo seemms to only be included in a given epoch when there are
% spikes during that epoch. 

% S_e3info = [e3info_unpacked{:}];

% Check that the number of cells are equal between the two arrays
if length(S_e3spikes) ~= length(e3info_unpacked)
    error("Length of arrays do not match")
end

for i = 1:length(S_e3spikes)
    S_e3spikes(i).tag = e3info_unpacked{i}.area;
end


figure;

subplot(2,2,1)
hold on;
title("Epoch 3 Raster All Cells")
for i=1:length(S_e3spikes)
    scatter(S_e3spikes(i).data(:,1), ones(size(S_e3spikes(i).data(:,1)))*i, ".")
end

subplot(2,2,2)
hold on;
title('Epoch 3 PFC Cells')
for i = 1:length(S_e3spikes)
    if S_e3spikes(i).tag == 'PFC'
        scatter(S_e3spikes(i).data(:,1), ones(size(S_e3spikes(i).data(:,1)))*i,".")
    end
end

subplot(2,2,3)
hold on;
title('Epoch 3 CA1 Cells')
for i = 1:length(S_e3spikes)
    if S_e3spikes(i).tag == 'CA1'
        scatter(S_e3spikes(i).data(:,1), ones(size(S_e3spikes(i).data(:,1)))*i,".")
    end
end



%% Looking at sleep states and LFP

% Decide what epoch and tetrode to look at. Must be odd numbered epoch for
% sleep.
epoch_ch = '11';
tetrode_ch = '25';
eeg_filetype = append('eeg01-',epoch_ch,'-',tetrode_ch);

load_rats = {'BG1_direct'}; % Which rats to load data from
filetype = {'ripples','sleep','rem','sws',eeg_filetype};

for i=1:length(load_rats)

    for k = 1:length(filetype)

        File_dir = dir(data_dir+"/"+load_rats(i)+"/**/*"+filetype(k)+"*");
        for j = 1:length(File_dir)
            load(string(fullfile(File_dir(j).folder, File_dir(j).name)));
        end

    end

end

% Extract proper epoch data
epoch_dbl = str2double(epoch_ch);
tetrode_dbl = str2double(tetrode_ch);
S_eeg = eeg{1,1}{1,epoch_dbl}{1,tetrode_dbl};
time_range = S_eeg.starttime:(1/S_eeg.samprate):S_eeg.endtime;
S_sleep = sleep{1,1}{1,epoch_dbl};
S_sws = sws{1,1}{1,epoch_dbl};
S_rem = rem{1,1}{1,epoch_dbl};
S_ripples = ripples{1,1}{1,epoch_dbl}{1,tetrode_dbl};


figure;
hold on;
title(sprintf("Sleep States | Epoch %d, Tet %d", epoch_dbl, tetrode_dbl))
plot(time_range, S_eeg.data)

sleep_labels = {};
for i = 1:length(S_sleep.starttime)
    x_vertices = [S_sleep.starttime(i),S_sleep.endtime(i),S_sleep.endtime(i),S_sleep.starttime(i)];
    y_vertices = [min(S_eeg.data),min(S_eeg.data),max(S_eeg.data),max(S_eeg.data)];
    patch(x_vertices, y_vertices, 'blue','FaceAlpha', 0.1, 'EdgeColor','none')

    if i == 1
        sleep_labels{i} = "sleep";
    else
        sleep_labels{i} = "";
    end
end

sws_labels = {};
for i = 1:length(S_sws.starttime)
    x_vertices = [S_sws.starttime(i),S_sws.endtime(i),S_sws.endtime(i),S_sws.starttime(i)];
    y_vertices = [min(S_eeg.data),min(S_eeg.data),max(S_eeg.data),max(S_eeg.data)];
    patch(x_vertices, y_vertices, 'green','FaceAlpha', 0.3, 'EdgeColor','none')

    if i == 1
        sws_labels{i} = "sws";
    else
        sws_labels{i} = "";
    end
end

rem_labels = {};
for i = 1:length(S_rem.starttime)
    x_vertices = [S_rem.starttime(i),S_rem.endtime(i),S_rem.endtime(i),S_rem.starttime(i)];
    y_vertices = [min(S_eeg.data),min(S_eeg.data),max(S_eeg.data),max(S_eeg.data)];
    patch(x_vertices, y_vertices, [0.8500 0.3250 0.0980],'FaceAlpha', 0.6, 'EdgeColor','none')

    if i == 1
        rem_labels{i} = "rem";
    else
        rem_labels{i} = "";
    end
end

ripple_labels = {};
for i = 1:length(S_ripples.starttime)
    x_vertices = [S_ripples.starttime(i),S_ripples.endtime(i),S_ripples.endtime(i),S_ripples.starttime(i)];
    y_vertices = [min(S_eeg.data),min(S_eeg.data),max(S_eeg.data),max(S_eeg.data)];
    patch(x_vertices, y_vertices, [0.4940 0.1840 0.5560],'FaceAlpha', 0.6, 'EdgeColor','none')

    if i == 1
        ripple_labels{i} = "ripples";
    else
        ripple_labels{i} = "";
    end
end

legend({"LFP", sleep_labels{:}, sws_labels{:}, rem_labels{:}, ripple_labels{:}})


%% Comparing LFP across tetrodes

% I want to look into ripples, but first I want to make sure that I am okay
% just referencing the LFP from a single tetrode. (Maybe I should just use
% the reference electrode for looking at LFP?)

load_rats = {'JS13_direct'}; % Which rats to load data from
filetype = {'eeg01-05'};

C_allfiles = {};

for i=1:length(load_rats)

        File_dir = dir(data_dir+"/"+load_rats(i)+"/**/*"+filetype(1)+"*");
        for j = 1:length(File_dir)
            S_file = load(string(fullfile(File_dir(j).folder, File_dir(j).name)));
            C_file = struct2cell(S_file);
            C_file = C_file{:}; % Remove extra cell outer layer
            C_allfiles(j) = C_file;
        end

end

% load information about the tetrodes


% So now looking at the LFP data on each tetrode for one epoch:


% figure;
% hold on;
% title("LFP on All Tetrodes for one Epoch")
% for i = 1:length(C_allfiles)
%     S_eeg = C_allfiles{1,i}{1,5}{1,i};
%     time_range = S_eeg.starttime:(1/S_eeg.samprate):S_eeg.endtime;
%     plot(time_range,S_eeg.data+1000*(i-1))
% 
% end


figure;
subplot(1,2,1)
hold on;
title("CA1 Tetrode LFP")
plotcount = 0;
for i = 1:length(C_allfiles)
    S_eeg = C_allfiles{1,i}{1,5}{1,i};
    if tetinfo{1,1}{1,5}{1,i}.area == 'CA1'
        time_range = S_eeg.starttime:(1/S_eeg.samprate):S_eeg.endtime;
        plot(time_range,S_eeg.data+1000*plotcount)
        plotcount = plotcount + 1;
    end
end

subplot(1,2,2)
hold on;
title("PFC Tetrode LFP")
plotcount = 0;
for i = 1:length(C_allfiles)
    S_eeg = C_allfiles{1,i}{1,5}{1,i};
    if tetinfo{1,1}{1,5}{1,i}.area == 'PFC'
        time_range = S_eeg.starttime:(1/S_eeg.samprate):S_eeg.endtime;
        plot(time_range,S_eeg.data+1000*plotcount)
        plotcount = plotcount + 1;
    end
end

%%

figure;
hold on;
S_eeg = C_allfiles{1,3}{1,5}{1,3};
time_range = S_eeg.starttime:(1/S_eeg.samprate):S_eeg.endtime;
plot(time_range,S_eeg.data)




