% dp: Douglas-Peucker, ma: moving average, cm: comb
% tolerance in meters (Douglas-Peucker)or minimum spacing between points
% for my moving average approach
function bedsout = sparsifyBedtraces(bedsin,tol,meth)

parser = inputParser;
addRequired(parser,'bedsin',@isnumeric)
addRequired(parser,'tol',@isscalar)
addRequired(parser,'meth',@ischar)

parse(parser,bedsin,tol,meth)
bedsin = parser.Results.bedsin;
tol = parser.Results.tol;
meth = parser.Results.meth;

% validate input
meth = validatestring(meth,{'cm','ma','dp'});
assert(size(bedsin,2)==4,'bed traces need 3 coordinates and one ID')

% parallelize, moving average is demanding
parll = 12;
par = gcp('nocreate');
if ~isempty(parll) && isempty(par)
    par = parpool(parll);
elseif ~isempty(parll) && par.NumWorkers ~= parll
    delete(par);
    par = parpool(parll);
elseif isempty(parll)
    delete(par);
end

% get all bed ids
beds = bedsin(:,4);
bedsunique = unique(beds);
nbeds = numel(bedsunique);
% get coordinates
pZ = bedsin(:,1:3);

% for each bed, perform some sort of point sparsification
% options are:
%   - Douglas-Peucker algorithm 'dp'
%   - moving average 'ma'
%   - combing 'cm' : this method uses the tolerance to find points that are
%   too close together, and it simply takes the first point in such a
%   cluster, as opposed to the moving average, which averages the point
%   locations for points that are within the tolerance range

% array of bedtraces after processing
switch meth
    case 'cm'
        parfor j = 1:nbeds
            disp(j)
            % points on current bed trace, will ideally dwindle
            curp = pZ(beds == bedsunique(j),:);
            % start at first indexed point
            curidx = 1;
            % start with 0 minimum spacing between points on bed trace
            minspace = 0;
            while minspace < tol
                % get all points within the tolerance distance of the
                % currently indexed point
                inidx = rangesearch(curp,curp(curidx,:),tol);
                % if there are points in the tolerance distance, take the
                % first
                if numel(inidx{1}) > 1
                    outidx = true(size(curp,1),1);
                    outidx(inidx{1}) = false;
                    % replace
                    % might be faster to switch this order?
                    curp = [curp(outidx,:); curp(curidx,:)];
                    % reset search index
                    curidx = 1;
                    % recompute minimum distance to reflect reduction in
                    % points
                    minspace = min(pdist(curp));
                else
                    curidx = curidx+1;
                    minspace = min(pdist(curp));
                end
            end
            bedsout{j} = [curp, bedsunique(j)*ones(size(curp,1),1)];
        end
    case 'dp'
        parfor j = 1:nbeds
            curp = pZ(beds == bedsunique(j),:);
            curps = dpsimplify(curp,tol);
            bedsout{j} = [curps, bedsunique(j)*ones(size(curps,1),1)];
        end
    case 'ma'
        parfor j = 1:nbeds
            disp(j)
            % points on current bed trace, will ideally dwindle
            curp = pZ(beds == bedsunique(j),:);
            % array to track how many points are combined into each
            % averaged point
            counts = ones(size(curp,1),1);
            % start at first indexed point
            curidx = 1;
            % start with 0 minimum spacing between points on bed trace
            minspace = 0;
            while minspace < tol
                % get all points within the tolerance distance of the
                % currently indexed point
                inidx = rangesearch(curp,curp(curidx,:),tol);
                if numel(inidx{1}) > 1
                    outidx = true(size(curp,1),1);
                    outidx(inidx{1}) = false;
                    % replace
                    % compute weighted sum of points in range
                    sumcounts = sum(counts(inidx{1}));
                    avgpt = sum( ...
                      bsxfun(@times,curp(inidx{1},:),counts(inidx{1}))...
                        ) / sumcounts;
                    % might be faster to switch this order?
                    curp = [curp(outidx,:); avgpt];
                    counts = [counts(outidx); sumcounts];
                    % reset search index
                    curidx = 1;
                    % recompute minimum distance to reflect reduction in
                    % points
                    minspace = min(pdist(curp));
                else
                    curidx = curidx+1;
                    minspace = min(pdist(curp));
                end
            end
%             bedsout = [bedsout; curp, bedsunique(j)*ones(size(curp,1),1)];
            bedsout{j} = [curp, bedsunique(j)*ones(size(curp,1),1)];
        end
end
bedsout = cat(1,bedsout{:});

end