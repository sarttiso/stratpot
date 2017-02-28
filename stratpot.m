% pZ: matrix of bedtraces with [x,y,z,bedID]
% bedID: vector of bed ids corresponding to coordinates in pZ
% pG: matrix of bedding attitude measurements with [x,y,z]
% G: matrix of gradient components [Gx,Gy,Gz]
% X: coordinates for which to krige stratigraphic potential [x,y,z]
% range: initial guess for range, or range(s) to be used. If you are
%   specifying the range(s) to be used, you must also specify
%   'usevarioparams' as true. For an initial guess, you may specify only
%   one value for range. If you specify only one value and then specify
%   'usevarioparams' as true, the value will be used for all gradient
%   components and directions. If you don't specify one value, you must
%   specify 1x3 or ndirx3 values. If you specify 1x3 ranges, then
%   the first is the range for Gx, second for Gy, third for Gz. If you
%   specfiy 1x3 ranges and have ndir > 1, then the same range will be use
%   in each direction for each gradient component. 
% sill: initial guess for sill(s). Operates as range above.
%
% OPTIONAL:
%
% parameters for variogram_emp, see description there. used if
% 'usevarioparams' is false (by default).
% 'ndir': (default 1, isotropic case)
% 'ndim': (default 0, isotropic case)
% 'nbins': (default 10)
% 'maxdist': (default empty)
% 'width': (default 30º)
% 'phi': (default 0)
%
% 'driftcor': (true/false) whether or not to correct data for drift 
%   (default false)
%
% 'driftfun': ('linear','quadratic') type of surface to fit for drift
%   correction of gradient data prior to selecting a variogram model;
%   different from drift of potential in kriging system
%
% 'nugget': initial guess for nugget(s). Operates like range and sill 
%   above. Not necessary to specify if 'usevarioparams' is true, since then 
%   only sill and range matter.
%
% 'usegui': true/false designating whether user can change
%   automatically fitted variogram parameters for each gradient component
%   and direction by using a GUI (default false). if true then
%   'usevarioparams' must be false.
%
% 'usevarioparams': true/false stipulating whether or not stratpot should 
%   use the user provided ranges, sills, and nuggets for each direction of
%   anisostropy (default false). If true, then 'adjustvarioparams' must be
%   false.
%
% 'anidir': directions of anisotropy, specified as azimuthal and
%   colatitudinal angles (like those in r as computed by variogram_emp).
%   if 'usevarioparams' is true, will by default be computed for 'ndir' and
%   'ndim' as in variogram_emp. if specified, must have ndir rows and
%   ndim-1 columns. For the isotropic case, 'anidir' is 1 (by default) 
%
% 'krigdrift': ('linear','quadratic'), default: 'linear'. Specifies the 
%   order of a polynomial basis for the kriging mean drift.
%
% 'potcovmod': default: 'cubic' specifies the model to use for potential
%   covariance.
%
% LIMITATIONS/TO DO:
%   - currently only considers 2-D anisotropy
%   - no support for faults yet
%   - add support for no drift in kriging matrix
%   - only implements dual form, so no consideration of kriging variance
%   - currently only supports cubic potential covariance model
%
% DEPENDENCIES:
%   - MultiPolyRegress(), by Ahmet Cecen
%   - Kriging toolbox, by Wolfgang Schwanghart and augmented by Adrian
%       Tasistro-Hart
%
% BUGS:
%   - no known bugs
%
% Adrian Tasistro-Hart Feb. 24 2017

function Zsampled = stratpot(pZ,bedID,pG,G,X,range,sill,varargin)


%% PARSING
parser = inputParser;

addRequired(parser,'pZ',@isnumeric);
addRequired(parser,'bedID',@isnumeric);
addRequired(parser,'pG',@isnumeric);
addRequired(parser,'G',@isnumeric);
addRequired(parser,'X',@isnumeric);
addRequired(parser,'range',@isnumeric);
addRequired(parser,'sill',@isnumeric);

% arguments for variogram_emp if 'usevarioparams' is false
addParameter(parser,'ndir',1,@isnumeric);
addParameter(parser,'ndim',0,@isnumeric);
addParameter(parser,'nbins',10,@isnumeric);
addParameter(parser,'maxdist',[],@isnumeric);
addParameter(parser,'width',30,@isnumeric);
addParameter(parser,'phi',0,@isnumeric);

addParameter(parser,'driftcor',false,@islogical);
% by default, if using drift correction, use planes
addParameter(parser,'driftfun','linear',@ischar);  
% nugget
addParameter(parser,'nugget',0,@isnumeric);
% whether or not user can use GUI to modify varioparams
addParameter(parser,'usegui',false,@islogical);
% force stratpot to use given ranges, sills, and nuggets. 
addParameter(parser,'usevarioparams',false,@islogical);
% specify directions of anisotropy for uservarioparams case
addParameter(parser,'anidir',[],@isnumeric);
% specify the order of kriging drift
addParameter(parser,'krigdrift','linear',@ischar);
% specify potential covariance model
addParameter(parser,'potcovmod','cubic',@ischar);

% parse here
parse(parser,pZ,bedID,pG,G,X,range,sill,varargin{:});

% required inputs
pZ =        parser.Results.pZ;
bedID =     parser.Results.bedID;
pG =        parser.Results.pG;
G =         parser.Results.G;
X =         parser.Results.X;
range0 =    parser.Results.range;
sill0 =     parser.Results.sill;

% variogram_emp parameters
ndir =      parser.Results.ndir;
ndim =      parser.Results.ndim;
nbins =     parser.Results.nbins;
maxdist =   parser.Results.maxdist;
wid =       parser.Results.width;
phi =       parser.Results.phi;

% strat pot parameters
driftcor =  parser.Results.driftcor;
driftfun =  parser.Results.driftfun;
nug0 =      parser.Results.nugget;
useb0 =     parser.Results.usevarioparams;
r =         parser.Results.anidir;
krigdrift = parser.Results.krigdrift;
potcovmod = parser.Results.potcovmod;
usegui =    parser.Results.usegui;

driftfun = validatestring(driftfun,{'linear','quadratic'});
krigdrift = validatestring(krigdrift,{'linear','quadratic'});
potcovmod = validatestring(potcovmod,{'cubic'});

% validate positions, gradients, bedIDs
assert(size(pG,2) == 3, 'pG must have x,y,z coordinates')
assert(size(pZ,2) == 3, 'pZ must have x,y,z coordinates')
assert(size(G,2) == 3, 'G must have x,y,z components')
assert(size(pG,1) == size(G,1), 'pG and G must have same number of rows')
assert(size(pZ,1) == size(bedID,1), ...
    'pZ and bedID must have same number of rows')

% check vars from argin
assert(((ndir == 1) && (ndim == 0)) || ...
       ((ndir > 1) && (ndim == 2 || ndim == 3)), ...
       'check ndir and ndmin');

% check given range(s), sill(s), nugget(s)
% range
% each gradient component has same range in each direction
if numel(range0) == 1  
    range0 = range0 * ones(ndir,3);
% range for each gradient component is same in every direction
elseif all(size(range0) == [1,3]) 
    range0 = repmat(range0,ndir,1);
else
    assert(all(size(range0) == [ndir,3]), ...
        'given range(s) must be 1x1, 1x3, or ndirx3');
end
% sill
if numel(sill0) == 1
    sill0 = sill0 * ones(ndir,3);
elseif all(size(sill0) == [1,3])
    sill0 = repmat(sill0,ndir,1);
else
    assert(all(size(sill0) == [ndir,3]), ...
        'given sill(s) must be 1x1, 1x3, or ndirx3');
end

%% PRE-PROCESSING

% demean position data 
xmean = mean(vertcat(pZ,pG,X));
pZ = bsxfun(@minus,pZ,xmean);
pG = bsxfun(@minus,pG,xmean);
X = bsxfun(@minus,X,xmean);

% make sure that gradient data is unit magnitude
G = bsxfun(@rdivide,G,sqrt(sum(G.^2,2)));

% input data size variables
m = size(G,1);
n = nincrements(bedID);
ntraces = length(bedID);

% bed variables
beds = unique(bedID);
nbeds = length(beds);

% not necessary to compute and fit empirical variograms if the user has 
% supplied necessary ranges, sills, and nuggets
if ~useb0
% check nugget
if numel(nug0) == 1
    nug0 = nug0 * ones(ndir,3);
elseif all(size(nug0) == [1,3])
    nug0 = repmat(nug0,ndir,1);
else
    assert(all(size(nug0) == [ndir,3]), ...
        'given nugget(s) must be 1x1, 1x3, or ndirx3');
end
% ncoorddist provides indices of columns to include when computing
% distances; for isotropic case or fully 3-D anisotropy, consider all
% coordinate axes when computing distances, but for 2-D anisotropy only
% consider the first two coordinate axes (x,y)
if ndim == 3 || ndim == 0
    ncoorddist = 3;
elseif ndim == 2
    ncoorddist = 2;
end
%% COMPUTE VARIOGRAMS, WITH OR WITHOUT DRIFT CORRECTION
% drift correction should depend on requested anisotropy. with 2-D
% anisotropy, we're just interested in drift with horizontal distance. with
% 3-D anisotropy (ndim = 3), or the isotropic case (ndim = 0), we must also
% consider drift vertically, since vertical displacement contributes to
% total displacement. 

    % DRIFT CORRECTION
    if driftcor
        % drift only horizontally
        if ndim == 2
            if strcmp(driftfun,'linear')
                drftx = MultiPolyRegress([pG(:,1),pG(:,2)],G(:,1),1);
                drfty = MultiPolyRegress([pG(:,1),pG(:,2)],G(:,2),1);
                drftz = MultiPolyRegress([pG(:,1),pG(:,2)],G(:,3),1);
            elseif strcmp(driftfun,'quadratic')
                drftx = MultiPolyRegress([pG(:,1),pG(:,2)],G(:,1),2);
                drfty = MultiPolyRegress([pG(:,1),pG(:,2)],G(:,2),2);
                drftz = MultiPolyRegress([pG(:,1),pG(:,2)],G(:,3),2);
            else
               error('driftfun not valid'); 
            end
        % drift horizontally and vertically
        elseif ndim == 3 || ndim == 0
            if strcmp(driftfun,'linear')
                drftx = MultiPolyRegress([pG(:,1),pG(:,2),pG(:,3)],...
                    G(:,1),1);
                drfty = MultiPolyRegress([pG(:,1),pG(:,2),pG(:,3)],...
                    G(:,2),1);
                drftz = MultiPolyRegress([pG(:,1),pG(:,2),pG(:,3)],...
                    G(:,3),1);
            elseif strcmp(driftfun,'quadratic')
                drftx = MultiPolyRegress([pG(:,1),pG(:,2),pG(:,3)],...
                    G(:,1),2);
                drfty = MultiPolyRegress([pG(:,1),pG(:,2),pG(:,3)],...
                    G(:,2),2);
                drftz = MultiPolyRegress([pG(:,1),pG(:,2),pG(:,3)],...
                    G(:,3),2);
            else
               error('driftfun not valid'); 
            end
        else
            error('check ndim for anisotropy')
        end
        % correct data with fitted surface
        Gxcor = drftx.Residuals;
        Gycor = drfty.Residuals;
        Gzcor = drftz.Residuals;

        [Vex,Hb,r,numobs] = variogram_emp(pG(:,1:ncoorddist),Gxcor,...
            'nbin',nbins,'ndir',ndir,'ndim',ndim,'maxdist',maxdist,...
            'width',wid,'phi',phi);
        Vey = variogram_emp(pG(:,1:ncoorddist),Gycor, ...
            'nbin',nbins,'ndir',ndir,'ndim',ndim,'maxdist',maxdist,...
            'width',wid,'phi',phi);
        Vez = variogram_emp(pG(:,1:ncoorddist),Gzcor, ...
            'nbin',nbins,'ndir',ndir,'ndim',ndim,'maxdist',maxdist,...
            'width',wid,'phi',phi);
    % NO DRIFT CORRECTION
    else
        [Vex,Hb,r,numobs] = variogram_emp(pG(:,1:ncoorddist),G(:,1),...
            'nbin',nbins,'ndir',ndir,'ndim',ndim,'maxdist',maxdist,...
            'width',wid,'phi',phi);
        Vey = variogram_emp(pG(:,1:ncoorddist),G(:,2), ...
            'nbin',nbins,'ndir',ndir,'ndim',ndim,'maxdist',maxdist,...
            'width',wid,'phi',phi);
        Vez = variogram_emp(pG(:,1:ncoorddist),G(:,3), ...
            'nbin',nbins,'ndir',ndir,'ndim',ndim,'maxdist',maxdist,...
            'width',wid,'phi',phi);
    end
    % concatenate gradient empirical variograms in 3-D array, 
    % nbins x ndir x 3
    V = cat(3,Vex,Vey,Vez);
    % ensure r includes both azimuthal and colatitudinal angles, might
    % want to add this processing to the end of variogram_emp (but
    % would then need to modify all related functions)
    if ndim == 2
        r = [r (pi/2)*ones(ndir,1)];
    end

%% FIT VARIOGRAM MODELS, get range and sill
% be very sure to remember whether or not you're working with drift
% corrected data.

% loop over all directions and fit sill and range to each direction for
% each gradient
sill_all = zeros(ndir,3);
range_all = zeros(ndir,3);
nugget_all = zeros(ndir,3);

% gui labels
gradlab = {'G_x','G_y','G_z'};
% need to update for 3-d anisotropy
if ndim == 0
    dirlab = {''};
elseif ndim == 2
    for j = 1:ndir
        dirlab{j} = sprintf('theta = %1.2f',r(j,1));
    end
else
    error('ndim must currenty be 0 or 2')
end

% gradient variogram functions
gfuns{1} = @gxx;
gfuns{2} = @gyy;
gfuns{3} = @gzz;

% consider each gradient component
for j = 1:3 
    % consider each direction
    for k = 1:ndir
        % isotropic
        if ndim == 0
            thetagfun = pi/4; % "average"
            phigfun = pi/2; % on average in plane, so keep z = 0
        % 2-D or 3-D anisotropy
        elseif ndim == 2 || ndim == 3
            thetagfun = r(k,1);
            phigfun = r(k,2);
        end
    
        fprintf('%s, %s ',gradlab{j},dirlab{k})
        % automatic (gradient descent) variogram fitting here
        [range_tmp, sill_tmp, nugget_tmp] = ...
            variogram_fit(Hb,V(:,k,j),range0(k,j),sill0(k,j),numobs,...
            'model',gfuns{j}, 'type','bounded','nugget',nug0(k,j),...
            'plotit',usegui);
        range_all(k,j) = range_tmp;
        sill_all(k,j) = sill_tmp;
        nugget_all(k,j) = nugget_tmp;
        fprintf('range = %1.3e, sill = %1.3e, nugget = %1.3e, C0 = %1.2f\n', ...
            range_all(k,j),sill_all(k,j),nugget_all(k,j),...
            sill_all(k,j)*range_all(k,j)^2/14);
        close(gcf)
    end
end

% if we use the given ranges and sills, impose that here
else
    range_all = range0;
    sill_all  = sill0;
    
    % need to check r is usevarioparams is true
    if isempty(r)
        if ndim == 0
            r = 1;
        elseif ndim == 2
            r = pointsOnCircle(ndir,0,180);
            r = [r(:,1), pi/2*ones(ndir,1)];
        elseif ndim == 3
            r = pointsOnSphere(ndir,0,[phi,90]);
            r = r(:,1:2);
        end
    else
        if ndir == 1
            assert(r==1,'for isotropic case, anidir is 1');
        else
            assert((size(r,1) == ndir) && (size(r,2) == ndim-1), ...
                'directions given for anidir are incorrect for given ndim and ndir')
        end
    end
end % end code for computing ranges, sills, and nuggets


%% CREATE COVARIANCE MATRICES
% use range(s) and sill(s) to fit models for potential variogram, consider 
% given anisotropy
% NEED TO ADD SUPPORT FOR NESTED MODELS
%       / KGxx      KGxy      KGxz      KGxZ-KGxZ'   F  \
%       | KGyx      KGyy      KGyz      KGyZ-KGyZ'      |
%  K =  | KGzx      KGzy      KGzz      KGzZ-KGzZ'      |
%       | KZGx-     KZGy-     KZGz-     KZZ - ...       |
%       |   KZGx'     KZGy'     KZGz'                   |
%
%

% potential covariance functions (cell array of functions)
cz   = potentialCovariance(range_all,sill_all,'model',potcovmod);
% remember that the derivative of the potential covariance function is
% alreay normalized by 1/r!!!
dcz  = potentialCovariance(range_all,sill_all,'derivative',1);
d2cz = potentialCovariance(range_all,sill_all,'derivative',2);

% initialize gradient covariance functions
cgxx = cgxxgen();
cgyy = cgyygen();
cgzz = cgzzgen();
cgxy = cgxygen();
cgxz = cgxzgen();
cgyz = cgyzgen();

% initialize gradient covariance matrices
Cgxx = zeros(m,m);
Cgyy = zeros(m,m);
Cgzz = zeros(m,m);
Cgxy = zeros(m,m);
Cgxz = zeros(m,m);
Cgyz = zeros(m,m);

% construct gradient covariance matrices from functions

% distance here does not depend on the ncoorddist variable because, as in
% Aug 2004-2005, the anisotropy model in the nested potential covariance
% function assumes that the axes of anisotropy are along the directions of
% the gradient components. Thus, there is the anisotropy that is captured
% in each of the ndir directions, and for each of these directions, there
% is a nested potential covariance model, and each of these nested
% potential covariance models assumes that the corresponding axes of
% anisotropy are along the coordinate axes.
for j = 1:m
    % gradients
    for k = 1:m
        dist = pG(j,:)-pG(k,:);
        % convert to spherical coordinates
        [theta,phi,~] = cart2sph2(dist);
        % based on inter-point angles choose best nested function by
        % finding closest direction of anisotropy
        aniidx = getAniidx(dist);
        
        Cgxx(j,k) = cgxx{aniidx}(theta,phi,dist);
        Cgyy(j,k) = cgyy{aniidx}(theta,phi,dist);
        Cgzz(j,k) = cgzz{aniidx}(theta,phi,dist);
        Cgxy(j,k) = cgxy{aniidx}(theta,phi,dist);
        Cgxz(j,k) = cgxz{aniidx}(theta,phi,dist);
        Cgyz(j,k) = cgyz{aniidx}(theta,phi,dist);
    end
end

% initialize potential-gradient cross-covariance functions
czgx = czgxgen();
czgy = czgygen();
czgz = czgzgen();

% initialize potential and potential-gradient cross-covariance matrices
Czgx = zeros(n,m);
Czgy = zeros(n,m);
Czgz = zeros(n,m);
Cz   = zeros(n,n);

% get indices for increments
secondIdx = zeros(ntraces,n);
firstIdx = zeros(ntraces,n);
for j = 1:n
    [secondIdx(:,j),firstIdx(:,j)] = incrementIdx(bedID,j,n);
end
secondIdx = logical(secondIdx);
firstIdx = logical(firstIdx);

% construct potential and potential-gradient covariance matrices
for j = 1:n
%     [rsecondIdx, rfirstIdx] = incrementIdx(bedID,j,n);
    % cross-covariance, gradient and increments 
    for k = 1:m
        % distance between gradient and jth point on bed trace
        dist2 = pZ(secondIdx(:,j),:) - pG(k,:);
        dist1 = pZ(firstIdx(:,j),:) - pG(k,:);
        % convert to spherical coordinates
        [theta2,phi2,~] = cart2sph2(dist2);
        [theta1,phi1,~] = cart2sph2(dist1);
        % get indices of nested covariance models to use
        aniidx2 = getAniidx(dist2);
        aniidx1 = getAniidx(dist1);
        Czgx(j,k) = czgx{aniidx2}(theta2,phi2,dist2) - ...
            czgx{aniidx1}(theta1,phi1,dist1);
        Czgy(j,k) = czgy{aniidx2}(theta2,phi2,dist2) - ...
            czgy{aniidx1}(theta1,phi1,dist1);
        Czgz(j,k) = czgz{aniidx2}(theta2,phi2,dist2) - ...
            czgz{aniidx1}(theta1,phi1,dist1);
    end  
    % increments
    for k = 1:n
        distrsjk = pZ(secondIdx(:,j),:) - pZ(secondIdx(:,k),:);
        distrs1k = pZ(secondIdx(:,j),:) - pZ(firstIdx(:,k),:);
        distrsj1 = pZ(firstIdx(:,j),:) - pZ(secondIdx(:,k),:);
        distrs11 = pZ(firstIdx(:,j),:) - pZ(firstIdx(:,k),:);
        Cz(j,k) = cz{getAniidx(distrsjk)}(distrsjk) - ...
            cz{getAniidx(distrs1k)}(distrs1k) - ...
            cz{getAniidx(distrsj1)}(distrsj1) + ...
            cz{getAniidx(distrs11)}(distrs11);
    end
    
end


%% DRIFT BASIS FUNCTIONS

if strcmp(krigdrift,'linear')
    % linear
    f{1} = @(x,y,z) x;
    f{2} = @(x,y,z) y;
    f{3} = @(x,y,z) z;
    dfx{1} = @(x,y,z) 1;
    dfx{2} = @(x,y,z) 0;
    dfx{3} = @(x,y,z) 0;
    dfy{1} = @(x,y,z) 0;
    dfy{2} = @(x,y,z) 1;
    dfy{3} = @(x,y,z) 0;
    dfz{1} = @(x,y,z) 0;
    dfz{2} = @(x,y,z) 0;
    dfz{3} = @(x,y,z) 1;
    ndrift = 3;
elseif strcmp(krigdrift,'quadratic')
    % quadratic
    f{1} = @(x,y,z) x;
    f{2} = @(x,y,z) y;
    f{3} = @(x,y,z) z;
    f{4} = @(x,y,z) x.^2;
    f{5} = @(x,y,z) y.^2;
    f{6} = @(x,y,z) z.^2;
    f{7} = @(x,y,z) x.*y;
    f{8} = @(x,y,z) x.*z;
    f{9} = @(x,y,z) y.*z;
    dfx{1} = @(x,y,z) 1;
    dfx{2} = @(x,y,z) 0;
    dfx{3} = @(x,y,z) 0;
    dfx{4} = @(x,y,z) 2*x;
    dfx{5} = @(x,y,z) 0;
    dfx{6} = @(x,y,z) 0;
    dfx{7} = @(x,y,z) y;
    dfx{8} = @(x,y,z) z;
    dfx{9} = @(x,y,z) 0;
    dfy{1} = @(x,y,z) 0;
    dfy{2} = @(x,y,z) 1;
    dfy{3} = @(x,y,z) 0;
    dfy{4} = @(x,y,z) 0;
    dfy{5} = @(x,y,z) 2*y;
    dfy{6} = @(x,y,z) 0;
    dfy{7} = @(x,y,z) x;
    dfy{8} = @(x,y,z) 0;
    dfy{9} = @(x,y,z) z;
    dfz{1} = @(x,y,z) 0;
    dfz{2} = @(x,y,z) 0;
    dfz{3} = @(x,y,z) 1;
    dfz{4} = @(x,y,z) 0;
    dfz{5} = @(x,y,z) 0;
    dfz{6} = @(x,y,z) 2*z;
    dfz{7} = @(x,y,z) 0;
    dfz{8} = @(x,y,z) x;
    dfz{9} = @(x,y,z) y;
    ndrift = 9;
end

%% DRIFT MATRICES

Fg = gradientDrift(pG,'order',krigdrift);
Fi = incrementDrift(pZ,bedID,'order',krigdrift);
F = [Fg, Fi];

% TO DO create fault functions?

%% CREATE KRIGING MATRIX 

CG = [Cgxx Cgxy' Cgxz'; Cgxy, Cgyy, Cgyz'; Cgxz, Cgyz, Cgzz];
CGI = [Czgx, Czgy, Czgz];
C = [CG, CGI'; CGI, Cz];
K = [C, F'; F, zeros(ndrift,ndrift)];

%% INVERT KRIGING MATRIX, GET WEIGHTS
B = [G(:,1); G(:,2); G(:,3); zeros(n,1); zeros(ndrift,1)];
w = K\B;

wx = w(1:m);                % weights on Czgx
wy = w(m+1:2*m);            % weights on Czgy
wz = w(2*m+1:3*m);          % weights on Czgz
wI = w(3*m+1:end-ndrift);   % weights on increments
wd = w(end-ndrift+1:end);   % weights on basis functions

%% INTERPOLATE AT REQUESTED POINTS

% solve kriging system at given points
Zsampled = krigeZ(X);

% save output


%% AUXILIARY FUNCTIONS

%% gradient component variograms
% variogram models for each gradient component, based on:
%   Chiles & Delfiner (1999) pg. 316
%   Aug & Chiles in Geostatistics (2004) pg 148-249
% pentaspheric & hole effect, first and second derivatives of the
% potential covariance function, taken to be cubic (see sources)
% need to redefine for each direction 
% alternatively, consider models presented in Lajaunie et al 1997
    % variogram of x component of gradient
    function gam = gxx(b,r)
        parsergxx = inputParser;
        addRequired(parsergxx,'b',@isnumeric);
        addRequired(parsergxx,'r',@isnumeric);

        parse(parsergxx,b,r);

        b = parsergxx.Results.b;
        r = parsergxx.Results.r;

        gam = b(2) * ...
        (sin(phigfun)^2*cos(thetagfun)^2 * ...
            (15/4*r./b(1)-  5*(r./b(1)).^3+9/4*(r./b(1)).^5) + ...
        (1-sin(phigfun)^2*cos(thetagfun)^2) * ...
            (15/8*r./b(1)-5/4*(r./b(1)).^3+3/8*(r./b(1)).^5));
    end
    % variogram of y component of gradient
    function gam = gyy(b,r)
        parsergyy = inputParser;
        addRequired(parsergyy,'b',@isnumeric);
        addRequired(parsergyy,'r',@isnumeric);

        parse(parsergyy,b,r);

        b = parsergyy.Results.b;
        r = parsergyy.Results.r;

        gam = b(2) * ...
        (sin(phigfun)^2*sin(thetagfun)^2 * ...
            (15/4*r./b(1)- 5*(r./b(1)).^3 +9/4*(r./b(1)).^5) + ...
        (1-sin(phigfun)^2*sin(thetagfun)^2) * ...
            (15/8*r./b(1)-5/4*(r./b(1)).^3+3/8*(r./b(1)).^5));
    end
    % variogram of z component of gradient
    function gam = gzz(b,r)
        parsergzz = inputParser;
        addRequired(parsergzz,'b',@isnumeric);
        addRequired(parsergzz,'r',@isnumeric);

        parse(parsergzz,b,r);

        b = parsergzz.Results.b;
        r = parsergzz.Results.r;

        gam = b(2) * ...
        (cos(phigfun)^2 * ...
            (15/4*r./b(1)- 5*(r./b(1)).^3 +9/4*(r./b(1)).^5) + ...
         sin(phigfun)^2 * ...
            (15/8*r./b(1)-5/4*(r./b(1)).^3+3/8*(r./b(1)).^5));
    end

%% kriging functions
    % krige potential
    function  Z = krigeZ(p)
        Z = zeros(size(p,1),1);
        for jj = 1:size(p,1)
            % five terms in interpolator
            terms = zeros(5,1);
            % gradient-increment cross-covariance
            for kk = 1:m
                distf = p(jj,:) - pG(kk,:);
                [thetaf,phif,~] = cart2sph2(distf);
                aniidxf = getAniidx(distf);
                terms(1) = terms(1) + ...
                    wx(kk)*czgx{aniidxf}(thetaf,phif,distf);
                terms(2) = terms(2) + ...
                    wy(kk)*czgy{aniidxf}(thetaf,phif,distf);
                terms(3) = terms(3) + ...
                    wz(kk)*czgz{aniidxf}(thetaf,phif,distf);
            end
            % increment covariance
            for kk = 1:n
                distf2 = p(jj,:)-pZ(secondIdx(:,kk),:);
                distf1 = p(jj,:)-pZ(firstIdx(:,kk),:);
                terms(4) = terms(4) + ...
                   wI(kk)* ( cz{getAniidx(distf2)}(distf2) - ...
                             cz{getAniidx(distf1)}(distf1) );
            end
            % drift
            for kk = 1:ndrift
                terms(5) = terms(5) + ...
                    wd(kk)*f{kk}(p(jj,1),p(jj,2),p(jj,3));
            end
            Z(jj) = sum(terms);
        end
    end

    % krige x component of gradient of potential
    function dZx = krigeDZx(p)
        dZx = zeros(size(p,1),1);
        for jj = 1:size(p,1)
            % five terms in interpolator
            terms = zeros(5,1);
            % gradient-increment cross-covariance
            for kk = 1:m
                distf = p(jj,:) - pG(kk,:);
                [thetaf,phif,~] = cart2sph2(distf);
                aniidxf = getAniidx(distf);
                terms(1) = terms(1) + ...
                    wx(kk)*cgxx{aniidxf}(thetaf,phif,distf);
                terms(2) = terms(2) + ...
                    wy(kk)*cgxy{aniidxf}(thetaf,phif,distf); 
                terms(3) = terms(3) + ...
                    wz(kk)*cgxz{aniidxf}(thetaf,phif,distf);
            end
            % increment covariance
            for kk = 1:n
                distf2 = p(jj,:)-pZ(secondIdx(:,kk),:);
                distf1 = p(jj,:)-pZ(firstIdx(:,kk),:);
                [thetaf2,phif2,~] = cart2sph2(distf2);
                [thetaf1,phif1,~] = cart2sph2(distf1);
                terms(4) = terms(4) + ...
                    wI(kk)*(czgx{getAniidx(distf2)}(thetaf2,phif2,distf2) - ...
                            czgx{getAniidx(distf1)}(thetaf1,phif1,distf1));
            end
            % drift
            for kk = 1:ndrift
                terms(5) = terms(5) + ...
                    wd(kk)*dfx{kk}(p(jj,1),p(jj,2),p(jj,3));
            end
            dZx(jj,1) = sum(terms);
        end
    end

    % krige y component of gradient of potential
    function dZy = krigeDZy(p)
        dZy = zeros(size(p,1),1);
        for jj = 1:size(p,1)
            % four terms in interpolator
            terms = zeros(5,1);
            % gradient-increment cross-covariance
            for kk = 1:m
                distf = p(jj,:) - pG(kk,:);
                [thetaf,phif,~] = cart2sph2(distf);
                aniidxf = getAniidx(distf);
                terms(1) = terms(1) + ...
                    wx(kk)*cgxy{aniidxf}(thetaf,phif,distf);
                terms(2) = terms(2) + ...
                    wy(kk)*cgyy{aniidxf}(thetaf,phif,distf);
                terms(3) = terms(3) + ...
                    wz(kk)*cgyz{aniidxf}(thetaf,phif,distf);
            end
            % increment covariance
            for kk = 1:n
                distf2 = p(jj,:)-pZ(secondIdx(:,kk),:);
                distf1 = p(jj,:)-pZ(firstIdx(:,kk),:);
                [thetaf2,phif2,~] = cart2sph2(distf2);
                [thetaf1,phif1,~] = cart2sph2(distf1);
                terms(4) = terms(4) + ...
                 wI(kk)*(czgy{getAniidx(distf2)}(thetaf2,phif2,distf2) - ...
                         czgy{getAniidx(distf1)}(thetaf1,phif1,distf1));
            end
            % drift
            for kk = 1:ndrift
                terms(5) = terms(5) + ...
                    wd(kk)*dfy{kk}(p(jj,1),p(jj,2),p(jj,3));
            end
            dZy(jj,1) = sum(terms);
        end
    end

    % krige z component of gradient of potential
    function dZz = krigeDZz(p)
        dZz = zeros(size(p,1),1);
        for jj = 1:size(p,1)
            % four terms in interpolator
            terms = zeros(5,1);
            % gradient-increment cross-covariance
            for kk = 1:m
                distf = abs(p(jj,:) - pG(kk,:));
                [thetaf,phif,~] = cart2sph2(distf);
                aniidxf = getAniidx(distf);
                terms(1) = terms(1) + ...
                    wx(kk)*cgxz{aniidxf}(thetaf,phif,distf);
                terms(2) = terms(2) + ...
                    wy(kk)*cgyz{aniidxf}(thetaf,phif,distf);
                terms(3) = terms(3) + ...
                    wz(kk)*cgzz{aniidxf}(thetaf,phif,distf);
            end
            % increment covariance
            for kk = 1:n
                distf2 = p(jj,:)-pZ(secondIdx(:,kk),:);
                distf1 = p(jj,:)-pZ(firstIdx(:,kk),:);
                [thetaf2,phif2,~] = cart2sph2(distf2);
                [thetaf1,phif1,~] = cart2sph2(distf1);
                terms(4) = terms(4) + ...
                 wI(kk)*(czgz{getAniidx(distf2)}(thetaf2,phif2,distf2) - ...
                         czgz{getAniidx(distf1)}(thetaf1,phif1,distf1));
            end
            % drift
            for kk = 1:ndrift
                terms(5) = terms(5) + ...
                    wd(kk)*dfz{kk}(p(jj,1),p(jj,2),p(jj,3));
            end
            dZz(jj,1) = sum(terms);
        end
    end

%% covariance functions

% gradient covariance functions, spherical coordinates, theta is angle
% counterclockwise from positive x-axis, phi is colatitude.
% DONT FORGET, first derivative of potential covariance is already
% normalized by 1/r
% h is the 3-components of displacement

    function cgxx = cgxxgen()
        cgxx = cell(ndir);
        for jj = 1:ndir
            cgxx{jj} = @(theta,phi,h) ...
                (sin(phi)^2*cos(theta)^2 - 1) * dcz{jj}(h) - ...
                (sin(phi)^2*cos(theta)^2) * d2cz{jj}(h);
        end
    end
    function cgyy = cgyygen()
        cgyy = cell(ndir);
        for jj = 1:ndir
            cgyy{jj} = @(theta,phi,h) ...
                (sin(phi)^2*sin(theta)^2 - 1) * dcz{jj}(h) - ...
                (sin(phi)^2*sin(theta)^2) * d2cz{jj}(h);
        end
    end
    function cgzz = cgzzgen()
        cgzz = cell(ndir);
        for jj = 1:ndir
            cgzz{jj} = @(theta,phi,h) ...
                (cos(phi)^2 - 1) * dcz{jj}(h) - ...
                (cos(phi)^2) * d2cz{jj}(h);
        end
    end
    function cgxy = cgxygen()
        cgxy = cell(ndir);
        for jj = 1:ndir
            cgxy{jj} = @(theta,phi,h) ...
                (sin(phi)^2*cos(theta)*sin(theta)) * ...
                (dcz{jj}(h) - d2cz{jj}(h));
        end
    end
    function cgxz = cgxzgen()
        cgxz = cell(ndir);
        for jj = 1:ndir
            cgxz{jj} = @(theta,phi,h) ...
                (sin(phi)*cos(theta)*cos(phi)) * ...
                (dcz{jj}(h) - d2cz{jj}(h));
        end
    end
    function cgyz = cgyzgen()
        cgyz = cell(ndir);
        for jj = 1:ndir
            cgyz{jj} = @(theta,phi,h) ...
                (sin(phi)*sin(theta)*cos(phi)) * ...
                (dcz{jj}(h) - d2cz{jj}(h));
        end
    end

% gradient-potential cross-covariance models
    function czgx = czgxgen()
        czgx = cell(ndir);
        for jj = 1:ndir
            czgx{jj} = @(theta,phi,h) ...
                -norm(h)*sin(phi)*cos(theta) * dcz{jj}(h);
        end
    end
    function czgy = czgygen()
        czgy = cell(ndir);
        for jj = 1:ndir
            czgy{jj} = @(theta,phi,h) ...
                -norm(h)*sin(phi)*sin(theta) * dcz{jj}(h);
        end
    end
    function czgz = czgzgen()
        czgz = cell(ndir);
        for jj = 1:ndir
            czgz{jj} = @(theta,phi,h) ...
                -norm(h)*cos(phi) * dcz{jj}(h);
        end
    end

%% find index of anisotropy when forming kriging matrices

    function aniidx = getAniidx(dist)
        if ndim == 0
            aniidx = 1;
        elseif ndim == 2 || ndim == 3
            % unit length vectors corresponding to each direction of
            % anisotropy
            us = sph2cart2([r,ones(ndir,1)]);
            % residuals
            res = zeros(ndir,1);
            for jj = 1:ndir
                u = us(jj,:)';
                res(jj) = norm(dist - u*dot(u,dist));
            end
            % get index of minimum residual
            [~,aniidx] = min(res);
        end
    end
end