% Script to analyze firing rate data of neurons in different states.
% Additionally this will be split by epoch to look for progressions over
% learning.
% Michael Satchell 11/02/23


%% Load data

clearvars;

dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
loadRats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01',
% 'spikes01','tetinfo','linfields01','rippletime01','pos01','nrninfo_MES','combstates_MES'
filetypes = {'nrninfo_MES','pos01','behavperform','combstates_MES'};

C_alldata = load_data(dataDir,loadRats,filetypes); % Load data to struct



if ~any(contains(filetypes,'pos01'))
    error("pos01 data must be loaded to run this analysis.")
end
if ~any(contains(filetypes,'behavperform'))
    error("behavperform data must be loaded to run this analysis.")
end
if ~any(contains(filetypes,'nrninfo_MES'))
    error("nrninfo_MES data must be loaded to run this analysis.")
end
if ~any(contains(filetypes,'combstates_MES'))
    error("combstates_MES data must be loaded to run this analysis.")
end

C_behper = C_alldata(find(contains(filetypes,'behavperform')),:);
C_nrninfo = C_alldata(find(contains(filetypes,'nrninfo_MES')),:);
C_pos = C_alldata(find(contains(filetypes,'pos01')),:);
C_combstates = C_alldata(find(contains(filetypes,'combstates_MES')),:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;
brainAreas = {'CA1','PFC'};




%% Split spikes up depending on the state they occur during

for r = 1:size(C_nrninfo,2)

    for e = 1:size(C_combstates{1,r},2)
        
        combOccs = C_combstates{1,r}{1,e}.combSortData;
        stateNames = C_comstates{1,r}{1,e}.stateNames;
        stateExist = ~cellfun(@isempty,C_combstates{1,4}{1,4}.sepData);

        
        for s = 1:size(stateNames,1)

            if stateExist(s) % If the state exists for this epoch

                

                for n = 1:size(C_nrninfo{1,r},1)
                    
                    S_nrn = C_nrninfo{1,r}{n,1};
        
        
                end
            end
        end

    end

end