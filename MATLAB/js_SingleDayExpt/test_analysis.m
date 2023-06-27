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

% The format for most of the data files is supposedly as follows:
% Day / Epoch(17) / Tetrode / Cell 

% There are 17 epochs, with 8 W track behavior sessions and 9 sleep
% sessions. Most animals have 32 tetrodes, although I think at least one
% has 64 tetrodes.


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
% spikes - Contains spike information for each cell. 

% More detailed information on the format of specific file types is as
% follows:
%
%
% linfields:
% day / epoch / tetrode / cell / trajectory
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
        filenames = {File_dir.name};
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

