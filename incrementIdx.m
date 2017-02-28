% this function returns the indices of the points associated with the jth
% increment of potential data for a given set of bed traces. the increments
% are always with respect to the first poing found on a bed trace.
% bedID: vector of IDs specifying which beds a set of points correspond to
% j: the increment under consideration
% nincrements: optional, specifies how many increments there are
% secondIdx: index of nonfixed point on bed trace
% firstIdx: index of fixed point on bed trace

function [secondIdx, firstIdx] = incrementIdx(bedID,j,varargin)
    parser = inputParser;
    addRequired(parser,'bedID',@isnumeric);
    addRequired(parser,'j',@isnumeric);
    addOptional(parser,'nincrements',0,@isnumeric);
    
    parse(parser,bedID,j,varargin{:});
    
    % bed IDs
    bedID = parser.Results.bedID;
    n = parser.Results.nincrements; 
    beds = unique(bedID);
    nbeds = length(beds);
    if n == 0
        n = nincrements(bedID);
    end
    
    m = length(bedID);      % number of points on bed traces
    
    % which increment to consider
    j = parser.Results.j;
    
    % if requested increment is too large, return error
    if j > n
        error('not sufficient increments for requested increment')
    end
    
    % get indices of first point on each bed trace
    firstIdx = zeros(m,1);
    for k = 1:nbeds
        firstIdx(find(bedID==beds(k),1)) = 1;
    end
    
    bedidx = 0;
    tmp = 0;
    while tmp < j
        bedidx = bedidx+1;
        rm = j - tmp;  % number of remaining indices to go until n
        % get total number of increments corresponding to current bed
        % (indexed by bedidx), keep looking until bed
        tmp = tmp + sum(ismember(bedID,beds(bedidx))) - 1;
    end
    clear tmp;
    
    % get index of fixed point on bed trace
    firstIdx = zeros(m,1);
    firstIdx(find(bedID==beds(bedidx),1)) = 1;
    
    % get indices of all other points on each bed trace
    secondIdx = ismember(bedID,beds(bedidx)) - firstIdx;
    
    % indices of nonzero entries in secondIdx
    idxs = find(secondIdx);
    
    % want index corresponding to leftover "space" before reaching n
    idx = idxs(rm);
    
    secondIdx = zeros(m,1);
    secondIdx(idx) = 1;
    
    firstIdx = logical(firstIdx);
    secondIdx = logical(secondIdx);
end