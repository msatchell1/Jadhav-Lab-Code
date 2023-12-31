function out = getspatialinfo(index, excludetimes, linpos, spikes, includestates, minV, varargin)
%   out = getspatialinfo(index, excludetimes, linpos, spikes, includestates,minV)
%          Returns the place information carried by the cell in bits using the
%          formula info = sum(P * (Ri/Rm) * log(Ri/Rm)) from Skaggs (1993)
%          The place information is computed based on the trajectory firing
%          rate maps
% position information = sum over all spatial bins of
%   (probability the animal was in that bin * (rate in that bin / average rate across bins) *
%   log_2(rate in that bin / average rate across bins)
% The probability the animal was in the bin is the fraction of total time
%   the animal was in the specific bin, and hopefully the rest is clear.


% add options
appendindex = 1;
peakrate = 0;
incltraj = [];
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'peakrate'
                peakrate = varargin{option+1};
            case 'incltraj'
                incltraj = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end


%caclulate trajdata if not specified
[state, lindist] = getbehavestate(linpos, index(1,1), index(1,2), includestates, 'minlinvelocity', minV);  %filters linpos.statematrix.traj by includestates
%find includetimes only
exclind = isExcluded(linpos{index(1,1)}{index(1,2)}.statematrix.time,excludetimes);  %find excluded time indeces
state(find(exclind)) = -1; %set excluded times to -1
%for each cell, compute trajdata
%cells = unique([ind(:,3:4); ind(:,5:6)],'rows');
trajdata = calclinfields(spikes,state,lindist,linpos, index);


%calculate spatial info for each traj
tempout = [];
if isempty(incltraj)
    incltraj = 1:length(trajdata);
end
for trj = incltraj
    if ~isempty(trajdata{trj}) && sum(isnan(trajdata{trj}(:,5))) < size(trajdata{trj},1)/2 & ...%if not empty & less than 1/2 bins are empty
        max(trajdata{trj}(:,5))>=peakrate %if traj has a peak on it
    %get probability animal was in each bin
        smocc = trajdata{trj}(:,2); %unsmoothed occupancy
        probocc = smocc./sum(smocc); %probability animal was in each bin, of all bins in traj
        
        %get rate in each bin
        normrate = trajdata{trj}(:,4)./nanmean(trajdata{trj}(:,4));  %unsmoothed occ normd firing rate / mean firing rate 
        
        if any(normrate)
            %get rid of bins with no firing
            incl = find(normrate);
            normrate = normrate(incl);
            probocc = probocc(incl);
            
            %get spatial info
            spinfoperbin = probocc .* normrate .* log2(normrate); %not summed over all bins
            
            %sum over all bins over all traj
            spinfo = nansum(spinfoperbin);
            
        else
            spinfo = NaN;
        end
%         if spinfo>20
%             keyboard
%         end
        
        tempout = [tempout; trj spinfo];
        
    end
end

%appendindex if desired
if appendindex
    out = [repmat(index, size(tempout,1), 1) tempout];
else
    out = tmpout;
end