function [infieldrate,outfieldrate]=mapinfo_get_inout_pfrates_2(spmap,dwmap,pfmap)
% get mean firing rates in and out of field
% get infield
% 31/05/2016 altered to work with threeD space program

pf_idx=find(pfmap > 0);
spsum=sum(spmap(pf_idx));
dwsum=sum(dwmap(pf_idx));
infieldrate=(spsum/dwsum); %*Fs;

% get outfield
pf_idx=find(pfmap==0);
spsum=sum(spmap(pf_idx));
dwsum=sum(dwmap(pf_idx));
outfieldrate=(spsum/dwsum); %*Fs;

