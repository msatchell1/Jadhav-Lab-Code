function [C_alldata_clipped] = clip_17_epochs(C_alldata)
%CLIP_17_EPOCHS Removes epoch numbers > 17 from input cell array.
% 
% Removes all data associated with epoch numbers above 17. This is because
% some rats have 19 epochs for some of the filetypes. Justin did not
% include this data in his analysis and so he said I should not either. 
% Data passed containing 17 or fewer epochs is unaltered.

% Inputs:
%   C_alldata: (num filetypes) x (num rats) cell array containing all loaded
%   data. Each cell is size 1 x (num epochs).
% 
% Outputs:
%   C_alldata_clipped: Same as C_alldata, but with the extra epochs
%   removed.
%
% Michael Satchell 07/11/23

C_alldata_clipped = {};

% Loop through C_alldata and copy to C_alldata_clipped, removing extra
% epochs. There might be more list-comprehension-like way to do this, but I
% don't know how to do it with a cell array.
for ft = 1:size(C_alldata,1)
   for r = 1:size(C_alldata,2)

       epoch_array = C_alldata{ft,r};

       if length(epoch_array) <= 17
        C_alldata_clipped{ft,r} = epoch_array;
       else
        C_alldata_clipped{ft,r} = epoch_array(1:17); % Selects first 17 epochs and stores them.
       end

   end
end


end

