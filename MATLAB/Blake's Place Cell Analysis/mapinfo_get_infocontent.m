function info=mapinfo_get_infocontent(spmap,dwmap,numrow,numcol)
% info=mapinfo_get_infocontent(spmap,dwmap,numrow,numcol);
info = 0;
dwsum=sum(dwmap(:));
if dwsum == 0; info=0; return; end
rate= sum(spmap(:))/dwsum;             % get firing rate again this time after usampling correction
if rate == 0; info=0; return; end

for i = 1:numrow
  for j = 1:numcol
    lambda = 0;
    if dwmap(i, j) > 0; lambda = (spmap(i, j) / dwmap(i, j)) / rate;end
    prob = dwmap(i, j) / dwsum;
    if not(lambda == 0); info = info + prob * lambda * log2(lambda);end
  end
end