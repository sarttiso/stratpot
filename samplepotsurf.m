% This function allows the user to randomly sample data patches along
% transects of miniumum length for data that exists on a 2-D surface in 2-D
% or 3-D and for which potential has been interpolated. The notion of a
% surface is important, since this function is appropriate only for
% datasets such as drone-derived point clouds of color, which sample a
% surface. A critical condition here is that the surface does
% not fold over on itself when seen from the z-axis. This is true of all
% natural surfaces of interest that a drone could image meaningfully.
% 
% In this situation, any surface can be "cookie-cut" along the z-axis to
% isolate a patch of a given width that does not fold back on itself. These
% patches have two endpoints in potential space the define the extent of
% the randomly sampled patch. The ultimate goal is then to form a user
% specified number of transects of user specified width in the original
% space (of drone imagery, for instance) and then return either all the
% data binned by potential for each transect or averages/variances/etc of
% the data along bins of potential. 
%
% These transects can then be used as input for a frequency domain
% averaging of time series, like pchave.
%
% Note that he input data need to be sampled from within a convex polygon,
% otherwise transects will pass through unsampled space.
%
% Note: 'patch' and 'transect' are used interchangeably here
% 
% IN:
% pS: Data locations in original space (ndata x ndim) ndim in [2,3]
% pZ: Interpolated potential for the data (ndata x 1)
% dat: Data values (ndata x nfields)
% npatch: Number of patches
% len: Minimum length of each patch in potential space. 
% wid: Width of each patch, in units of the coordinates given for the data
%   locations.
% dz: Sample spacing in potential space. If it is too small (i.e. no data
%   falls within bins of potential for some patches) the user will be
%   warned.
% 'output': Specify what the output should be for each patch.
%       'average' (default): average the input data within the bins of
%           potential over each patch.
%       'std': return the standard deviation of input data within bins of
%           potential in each patch
%       'bins': return the indices of the original data into bins, along
%           with the bin edges for each patch.
%
% OUT:
% Y: Cell array of length npatch with each element containing an array
%   corresponding to a patch and containing the requested output
% T: Cell array of length npatch with each element containing an array of
%   coordinates (or bin edges for 'output'->'bins') in potential space 
%   corresponding to the output data in Y.
% ptch: Cell array containing the patch boundary coordinates in the
%   original space.
%
% REQUIRES:
% - geom2D (David Legland's library)
%
% TO DO:
% - add ability to apply arbitrary number of functions to data in bins,
%   rather than having to select either avg or std.
% - might be good to fix seed to achieve deterministic output (i.e. on
%   randperm, rand)
% -add the ability to specify deterministic transect lengths rather than a
%   minimum length
%
% Adrian Tasitro-Hart, March 10 2018, adrianraph-at-gmail.com

function [Y,T,ptch] = samplepotsurf(pS,pZ,dat,npatch,len,wid,dz,varargin)

%% parse and validate

parser = inputParser;
addRequired(parser,'pS',@isnumeric)
addRequired(parser,'pZ',@isnumeric)
addRequired(parser,'dat',@isnumeric)
addRequired(parser,'npatch',@isscalar)
addRequired(parser,'len',@isscalar)
addRequired(parser,'wid',@isscalar)
addRequired(parser,'dz',@isscalar)
addParameter(parser,'output','average',@ischar)

% parse inputs
parse(parser,pS,pZ,dat,npatch,len,wid,dz,varargin{:});
pS     = parser.Results.pS;
pZ     = parser.Results.pZ;
dat    = parser.Results.dat;
npatch = parser.Results.npatch;
len    = parser.Results.len;
wid    = parser.Results.wid;
dz     = parser.Results.dz;
output = parser.Results.output;

% make sure everything is the right size
ndata = size(pS,1);     % number of data points
ndims = size(pS,2);     % dimensionality of data coordinates
nflds = size(dat,2);    % number of data fields (color, slope, etc.)
assert(ismember(ndims,[2,3]),'data must be 2-D or 3-D')
assert(length(pZ) == ndata,'pZ must be same length as pS')
assert(ndata == size(dat,1),'dat must be same length as pS')
% make sure pZ is column
pZ = pZ(:);

% make sure that requested length isn't too big for the span of the input
% data in potential space
assert(len < range(pZ),'len is too long for the given input data')

% validate output type
valid_output = {'average','std','bins'};
output = validatestring(output,valid_output);

%% generate patch endpoints

% 1) first find eligible endpoints for the selected length
% first randomize the potential coordinates to get a random sampling
randidx = randperm(ndata);
% now find the first endpoints by taking the first nptch eligible points,
% i.e. points that have other points further than len away from them
pt1idx = zeros(npatch,1);
c1 = 1; % counter 1, successful point 1 indices
c2 = 1; % counter 2, iterates through randidx
while c1 <= npatch
    if any(abs(pZ(randidx(c2))-pZ) > len)
        pt1idx(c1) = randidx(c2);
        c1 = c1+1;
    end
    c2 = c2+1;
end
clear c1 c2

% 2) find the second set of eligible endpoints, one for each of the
% previously found points
pt2idx = zeros(npatch,1);
% we will use this to generate random indices into the points that are far
% enough away to be an endpoint
tmprnd = rand(npatch,1);
parfor ii = 1:npatch
    tmpidx = find(abs(pZ(pt1idx(ii))-pZ) > len);
    pt2idx(ii) = tmpidx(ceil(tmprnd(ii)*length(tmpidx)));
end
clear tmprnd

% ensure that pt1 < pt2 in potential coordinates
tmpidx = pZ(pt1idx) > pZ(pt2idx); % find points for which this isnt true
tmppt = pt2idx(tmpidx);
pt2idx(tmpidx) = pt1idx(tmpidx);
pt1idx(tmpidx) = tmppt;
clear tmppt tmpidx

%% create patch slices and gather data into patches
% patch coordinates will exist in the original space of the data
% **uses geom2D, David Legland's library

% 1) first create patches
% cell array with each entry containing the four corners of every slice
ptch = cell(npatch,1);
parfor ii = 1:npatch
    curpt1 = pS(pt1idx(ii),1:2);
    curpt2 = pS(pt2idx(ii),1:2);
    mainl = createLine(curpt1,curpt2);
    orthl1 = orthogonalLine(mainl,curpt1);
    orthl2 = orthogonalLine(mainl,curpt2);
    corns1 = intersectLineCircle(orthl1,[curpt1(1),curpt1(2),wid/2]);
    corns2 = intersectLineCircle(orthl2,[curpt2(1),curpt2(2),wid/2]);
    corns = [corns1;corns2];
    % need to organize corners so that they enclose the patch, use convhull
    K = convhull(corns);
    ptch{ii} = corns(K(1:end-1),:);
end

% 2) now gather the data for the appropriate output
edges = cell(npatch,1);
T = cell(npatch,1);
inpatchidx = cell(npatch,1);    % indices of points in a given patch
binsidx = cell(npatch,1);       % indices of bins for points in a patch
emptbinidx = cell(npatch,1);    % indices of empty bins in a patch
parfor ii = 1:npatch
    % first, get all points in current patch
    curpatch = ptch{ii};
    inpatchidx{ii} = inpolygon(pS(:,1),pS(:,2),curpatch(:,1),curpatch(:,2));
    % now compute edges
    edges{ii} = min(pZ(inpatchidx{ii})):dz:max(pZ(inpatchidx{ii})); % careful with max here
    edges{ii}(end) = max(pZ(inpatchidx{ii}));
    % compute T on bin centers
    T{ii} = edges{ii}(1:end-1)+dz/2;
    % indices into bins relative to the data in the patch
    binsidx{ii} = discretize(pZ(inpatchidx{ii}),edges{ii});
    % check if any bins are empty
    emptbinidx{ii} = accumarray(binsidx{ii},1) == 0;
end
% if any bins are empty, sum over the total number of empty bins versus
% total bins and tell the user.
emptbins = 0;
totbins = 0;
for ii = 1:npatch
    emptbins = sum(emptbinidx{ii}) + emptbins;
    totbins = length(T{ii}) + totbins;
end
if emptbins > 0
    fprintf('%d empty bins of %d total bins\n',emptbins,totbins)
    % warn the user if more than 5% of bins are empty
    if emptbins/totbins > 0.05
        warning('more than 5% of bins are empty!')
    end
end

%% operate on the data in the patches
if ~strcmp(output,'bins')
    if strcmp(output,'average')
        fun = @mean;
    elseif strcmp(output,'std')
        fun = @std;
    end
    Y = cell(npatch,1);
    for ii = 1:npatch
        % for each data field, compute the desired output
        for jj = 1:nflds
            Y{ii}(:,jj) = accumarray(binsidx{ii},...
                                     dat(inpatchidx{ii},jj),[],fun);
        end
        % need to interpolate over bins that have no data here
        curidx = emptbinidx{ii};
        if sum(curidx > 0)
            % interpolate every data field
            for jj = 1:nflds
                Y{ii}(curidx,jj) = interp1(T{ii}(~curidx), ...
                    Y{ii}(~curidx,jj),T{ii}(curidx),'pchip');
            end
        end
    end
else
    % if user wants bin indices, need to recompute them relative to all of
    % the input data
    for ii = 1:npatch
        tmpidx = NaN(ndata,1);
        tmpidx(inpatchidx{ii}) = binsidx{ii};
        Y{ii} = tmpidx;
    end
end

end