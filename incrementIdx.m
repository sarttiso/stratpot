% this function returns the indices of points on bed traces that form
% increments for the interpolation of potential with stratpot. by default,
% the function returns two vectors (of equal lengths) of linear indices,
% each corresponding entry in these vectors indexing two points that form 
% a null increment of potential along a bed trace. the first vector
% contains the indices corresponding to points along a bed trace that
% change such that the increments along bed traces are defined by a point
% that is held constant on each bed trace. the indices of the points that
% are held constant on each bed are in the second vector. this approach
% corresponds to the case of taking x0 for defining the increment in Chiles
% 2004 (pg 314). if only the indices to one increment are desired, then
% these can also be requested.

% IN:
% bedID: vector of IDs specifying which beds a set of points correspond to
% nincrements: optional, specifies how many increments there are
% 'increment': (default empty) parameter specifying which increment to
%   return. By default, this is empty, and so the indices to all increments
%   are return. If the user desires one particular increment, then the
%   'increment',j jth increment will be returned, so long as j is less than
%   the total number of increments.
% OUT:
% secondIdx: linear indices of points relative to a fixed point on a bed 
%   trace, which define increments. If only one increment is requested,
%   then secondIdx is a vector of logical value with one true value
%   selecting the index of the required point for the requested increment.
% firstIdx: same as secondIdx, except that this vector indexes points that
%   are held constant on each bed trace, forming an anchor against which
%   the increment is defined.
% 
% Adrian Tasistro-Hart Feb 28 2017

function [secondIdx, firstIdx] = incrementIdx(bedID,varargin)
    parser = inputParser;
    addRequired(parser,'bedID',@isnumeric);
    addOptional(parser,'nincrements',0,@isnumeric);
    addParameter(parser,'increment',[],@isnumeric);
    
    parse(parser,bedID,varargin{:});
    
    % bed IDs
    bedID = parser.Results.bedID;
    n = parser.Results.nincrements; 
    beds = unique(bedID);
    nbeds = length(beds);
    if n == 0
        n = nincrements(bedID);
    end
    
    m = length(bedID);      % number of points on bed traces
    
    % which increment to consider if any
    inc = parser.Results.increment;
    
    % if user wants all increment indices
    if isempty(inc)
        % get indices for each increment
        secondIdx = zeros(n,1);
        firstIdx = zeros(n,1);
        for j = 1:n
            [secondTmp,firstTmp] = jthincrement(j);
            secondIdx(j) = find(secondTmp);
            firstIdx(j) = find(firstTmp);
        end
        % get row of index in each column
    % if user has requested specific increment, return it
    else 
        % if requested increment is too large, return error
        if inc > n
            error('not sufficient increments for requested increment')
        end
        [secondIdx, firstIdx] = jthincrement(inc);
    end
    
    
    %% AUX FUNCTION
    function [secondIdx, firstIdx] = jthincrement(j)
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
end