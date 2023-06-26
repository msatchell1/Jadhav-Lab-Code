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
% The eegref files are the reference LFP files. The delta, gamma, ripple,
% and theta files are filtered versions of the eeg files. Justin said that
% only the ripple data files need to have the reference LFP trace
% subtracted from it because noise is usually only high at that frequency.
%
% Ripples are calculated from tetrodes with cells, and are usually coherent
% across all tetrodes so it doesn't really matter which tetrode you use to
% measure the ripples. Ripples are also only detected for low velocity
% periods.


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


%% Load data

data_dir = '/mnt/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
load_rats = {'BG1_direct'}; % Which rats to load data from
filetype = 'ripple01-03-30'; % File type to load for each rat. All files will 
% be loaded that contain filetype in the file name. Use 'eeg01' when you
% want all the LFP files, otherwise eegref will also be included.



for i=1:length(load_rats)
    
    File_dir = dir(data_dir+"/"+load_rats(i)+"/**/*"+filetype+"*");
    filenames = {File_dir.name};
    for j = 1:length(File_dir)
        load(string(fullfile(File_dir(j).folder, File_dir(j).name)));
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

S_eeg = eeg{1,1}{1,3}{1,30};
time_range = eeg_struct.starttime:(1/eeg_struct.samprate):eeg_struct.endtime;

figure;
hold on;
title("Raw LFP Trace")
ylabel("Voltage (mV)??")
xlabel("Time (s)")
plot(time_range,eeg_struct.data)


