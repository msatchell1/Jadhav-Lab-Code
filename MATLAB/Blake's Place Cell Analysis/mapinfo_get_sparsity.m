function sparsity=mapinfo_get_sparsity(spmap,dwmap,numrow,numcol)
% calculate sparsity measure as per Skaggs and McNaughton (1996)

% get mean rate 
dwsum=sum(sum(dwmap));
spsum=sum(sum(spmap));
if dwsum > 0 ; rate=(spsum/dwsum);end;

%now process each pixel
sum1 = 0;
sum2 = 0;
for i = 1:numrow 
  for j = 1:numcol
    if dwsum > 0;prob = dwmap(i, j) / dwsum;end
    if dwmap(i, j) > 0; pixrate = spmap(i, j) / dwmap(i, j); end
    if dwmap(i, j) > 0; sum1 = sum1 + prob * pixrate; end
    if dwmap(i, j) > 0; sum2 = sum2 + prob * pixrate * pixrate; end
  end %j
end %i
if sum2 > 0 
  sparsity = (sum1 * sum1) / sum2;
else
 sparsity = 0;
end

