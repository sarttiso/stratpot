% given a set of bed traces at points p with corresponding bedIDs, generate
% the drift matrix of a given order for the potential increments
% p: npositions (rows) x ncoordinates (columns)
% bedID: npositions (rows) x 1, identifies bed for each of the given points
% 'order' : 'linear' or 'quadratic', defaults to 'linear'

function Fi = incrementDrift(p,bedID,varargin)
    parser = inputParser;
    addRequired(parser,'p',@isnumeric);
    addRequired(parser,'bedID',@isnumeric);
    addParameter(parser,'order','linear',@ischar);
    
    parse(parser,p,bedID,varargin{:});
    
    % bed trace data locations 
    p = parser.Results.p;
    m = size(p,1);  % number of input data
    D = size(p,2);  % dimension of input data space
    
    % bed IDs
    bedID = parser.Results.bedID;
    beds  = unique(bedID);
    nbeds = length(beds);   % total number of beds
    nk = hist(bedID,beds)'; % number of points for each bed
    ncomb = nk-1;           % number of point combinations per bed
    % check that all bed traces have at least two points
    if any(nk <= 1)
        error('bed with only one point on it, must have at least two'); 
    end
    n = sum(ncomb);  % total number of interpoint combinations
    
    % drift function order 
    order = parser.Results.order;
    order = validatestring(order,{'linear','quadratic'}); 
    
    % create drift matrix
    ncomb = cumsum([0; ncomb]); % easier to index regions of drift matrix
    if strcmp(order,'linear')
        Fi = zeros(D,n);
        for j = 1:n
            [secondIdx, firstIdx] = incrementIdx(bedID,n,'increment',j);
            for k = 1:D
                Fi(k,j) = p(secondIdx,k)-p(firstIdx,k);
            end
        end
        
    elseif strcmp(order,'quadratic') && D == 2
        Fi = zeros(5,n);
        for j = 1:n
            [secondIdx, firstIdx] = incrementIdx(bedID,j);
            % fill rows 1-4
            for k = 1:2
                Fi(k,j)   = p(secondIdx,k)-p(firstIdx,k);
                Fi(k+2,j) = p(secondIdx,k)^2-p(firstIdx,k)^2;
            end
            % fill 5th row
            Fi(5,j) = ...
             p(secondIdx,1)*p(secondIdx,2) - p(firstIdx,1)*p(firstIdx,2);
        end
        
    elseif strcmp(order,'quadratic') && D == 3
        Fi = zeros(9,n);
        for j = 1:n
            [secondIdx, firstIdx] = incrementIdx(bedID,j);
            % fill rows 1-6
            for k = 1:3
                Fi(k,j) = p(secondIdx,k)-p(firstIdx,k);
                Fi(k+3,j) = p(secondIdx,k)^2-p(firstIdx,k)^2;
            end
            % fill rows 7-9
            Fi(7,j) = ...
               p(secondIdx,1)*p(secondIdx,2) - p(firstIdx,1)*p(firstIdx,2);
            Fi(8,j) = ...
               p(secondIdx,1)*p(secondIdx,3) - p(firstIdx,1)*p(firstIdx,3);
            Fi(9,j) = ...
               p(secondIdx,2)*p(secondIdx,3) - p(firstIdx,2)*p(firstIdx,3);
        end
    else
        error('check input data')
    end
end