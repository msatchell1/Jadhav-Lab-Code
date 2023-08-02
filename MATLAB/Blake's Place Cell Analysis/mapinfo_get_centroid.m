function [cx,cy]=mapinfo_get_centroid(frmap,numrow,numcol)
% get centroid of place field map
xsum=sum(frmap);
ysum=sum(frmap');
xsum_weighted= xsum .* (1:numcol);
ysum_weighted= ysum .* (1:numrow);

cx=xsum_weighted/xsum;
cy=ysum_weighted/ysum;