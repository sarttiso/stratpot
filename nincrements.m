% given input bedIDs, returns the number of increments per bed

function [n, ncomb] = nincrements(bedID)
    
    beds  = unique(bedID);
    nk = hist(bedID,beds)'; % number of points for each bed
    
    % check that all bed traces have at least two points
    if any(nk <= 1)
        error('bed with only one point on it, must have at least two'); 
    end
    
    ncomb = nk-1;           % number of combinations per bed
    n = sum(ncomb);
    
end