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
        
        % combOccs = C_combstates{1,r}{1,e}.combSortData;
        stateNames = C_combstates{1,r}{1,e}.stateNames;

        % The index number into stateNames of which states exist for this
        % epoch.
        stateIdxs = find(~cellfun(@isempty,C_combstates{1,r}{1,e}.sepData));

        tl = tiledlayout(size(stateIdxs,1),size(stateIdxs,1));
        title(tl, sprintf("%s Epoch %i",loadRats{r},e))

        for s1 = stateIdxs'
            % Loop through a second set of states to get every combination
            for s2 = stateIdxs'

                s1Occs = C_combstates{1,r}{1,e}.sepData{s1,1};
                s2Occs = C_combstates{1,r}{1,e}.sepData{s2,1};
                s1s2Occs = [s1Occs; s2Occs]; % Combine so I can loop easier later.

                nexttile
                hold on
                
                xlabel(sprintf("FR in %s (Hz)",stateNames{s1}))
                ylabel(sprintf("FR in %s (Hz)",stateNames{s2}))
                set(gca, 'XScale', 'log', 'YScale', 'log');
                plot([0.01,40],[0.01,40],'k')
                % It's important to note that neurons in nrninfo_MES
                % without an area field also have no spike data in any
                % epochs. There are actually quite a few cells like this in
                % every rat.
                for n = 1:size(C_nrninfo{1,r},1)
                    
                    S_nrn = C_nrninfo{1,r}{n,1};
                    if S_nrn.eHasSpikeData(1,e) % If nrn has spike data for this epoch 
                        % if strcmp(S_nrn.type,"Int")
                        if strcmp(S_nrn.area,"PFC")
                            nrnColor = [152,8,230]/230; %"red";
                        % elseif strcmp(S_nrn.type,"Pyr")
                        elseif strcmp(S_nrn.area,"CA1")
                            nrnColor = [8,52,230]/230; %"blue";
                        end

                        spikes = S_nrn.eSpikeData{1,e}(:,1);
                        % Holds the number of spikes for a given neuron
                        % that occur during each state occurrance.
                        spikesInOcc = zeros(size(s1s2Occs,1),1);

                        % Loop through occurances of the two states s1 and
                        % s2, counting the number of spikes that fall into
                        % each occurrance during this epoch.
                        for o = 1:size(s1s2Occs,1)
                            isInOcc = (s1s2Occs(o,1) <= spikes) & (s1s2Occs(o,2) > spikes);
                            spikesInOcc(o,1) = sum(isInOcc);
                        end

                        % Sum respective spikes for each state
                        s1Sum = sum(spikesInOcc(s1s2Occs(:,3)==s1));
                        s2Sum = sum(spikesInOcc(s1s2Occs(:,3)==s2));
                        % total durations of states
                        s1Dur = sum(s1Occs(:,2) - s1Occs(:,1));
                        s2Dur = sum(s2Occs(:,2) - s2Occs(:,1));

                        scatter(s1Sum/s1Dur,s2Sum/s2Dur,MarkerFaceColor=nrnColor, ...
                            MarkerEdgeColor=nrnColor, Marker='.')
                        

                    end
                end
            end
        end
    pause
    end

end