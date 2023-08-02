function zcohere=mapinfo_get_spatial_coherence2(spmap,dwmap,numrow,numcol)
% calculate spatial coherence measure based on Kubie et al (1990)JN 10:1110
% calculate pair of values for each pixel. Ignore nonsampled pixels
% modified from mapinfo_get_spatial_coherence to cope with zero dwell

% actually works with 0s now - Blake DOUBLE SUM THE DWMAP 

if sum(sum(dwmap)) > 0
  numpairs=0;
  idx=find(dwmap==0);
  spmap(idx)=0;           % first ensure no spikes appear where dwell is 0
  for i = 2:numrow-1
    for j = 2:numcol-1
      %calculate firing rate for center pixel
      if dwmap(i, j) > 0
        pixrate = spmap(i, j) / dwmap(i, j);
        numpairs = numpairs + 1;
        coherarray(numpairs, 1) = pixrate;
        %calculate mean firing rate for neighbours
        matesp = spmap(i - 1, j - 1) + spmap(i, j - 1) + spmap(i + 1, j - 1) + spmap(i - 1, j) + spmap(i + 1, j) + spmap(i - 1, j + 1) + spmap(i, j + 1) + spmap(i + 1, j + 1);
        matedw = dwmap(i - 1, j - 1) + dwmap(i, j - 1) + dwmap(i + 1, j - 1) + dwmap(i - 1, j) + dwmap(i + 1, j) + dwmap(i - 1, j + 1) + dwmap(i, j + 1) + dwmap(i + 1, j + 1);
        if matedw > 0
          materate = matesp / matedw;
          coherarray(numpairs, 2) = materate; % this was at end before
        else
          numpairs = numpairs-1;     % not sampled therefore don't use
        end
        
      end % if
    end   % j
  end     % i
%now correlate pixel fr with that of neighbors 

  cohere = corr(coherarray(:,1),coherarray(:,2));
  if and(cohere > -0.999, cohere  < 0.999) 
    zcohere = 0.5 * log((1 + cohere) / (1 - cohere));      % z transformed coherence
  else
    zcohere=-9998;
  end
else
  zcohere=-9999;
end
