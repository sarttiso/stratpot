function gradout = sparsifyGradients(gradin,tol,meth)

parser = inputParser;
addRequired(parser,'gradin',@isnumeric)
addRequired(parser,'tol',@isscalar)
addRequired(parser,'meth',@ischar)

parse(parser,gradin,tol,meth)
gradin = parser.Results.gradin;
tol = parser.Results.tol;
meth = parser.Results.meth;

% validate input
meth = validatestring(meth,{'cm','ma','dp'});
assert(size(gradin,2)==6,'bed traces need 3 coordinates and 3 components')

% get coordinates
pG = gradin(:,1:3);
% get gradient values
G = gradin(:,4:6);
m = size(pG,1);

% perform some sort of point sparsification
% options are:
%   - moving average (not done yet)

% array of bedtraces after processing
switch meth
    case 'ma'
        % initial gradient locations, will ideally dwindle
        curp = pG;
        % initial gradient values
        curG = G;
        % array to track how many points are combined into each
        % averaged point
        counts = ones(m,1);
        % start at first indexed point
        curidx = 1;
        % start with 0 minimum spacing between points on bed trace
        minspace = 0;
        while minspace < tol
            % get all points within the tolerance distance of the
            % currently indexed point
            inidx = rangesearch(curp,curp(curidx,:),tol);
            % if we find points in range, average them
            if numel(inidx{1}) > 1
                outidx = true(size(curp,1),1);
                outidx(inidx{1}) = false;
                % replace
                % compute weighted sum of points in range
                sumcounts = sum(counts(inidx{1}));
                avgpt = sum( ...
                  bsxfun(@times,curp(inidx{1},:),counts(inidx{1}))...
                    ) / sumcounts;
                avgG = sum( ...
                    bsxfun(@times,curG(inidx{1},:),counts(inidx{1})) ...
                    ) / sumcounts;
                curp = [curp(outidx,:); avgpt];
                curG = [curG(outidx,:); avgG];
                counts = [counts(outidx); sumcounts];
                % reset search index
                curidx = 1;
                % recompute minimum distance to reflect reduction in
                % points
                minspace = min(pdist(curp));
            % otherwise, look at next point
            else
                curidx = curidx+1;
            end
        end
        % make sure gradients have unit length
        curG = bsxfun(@rdivide,curG,sqrt(sum(curG.^2,2)));
        gradout = [curp, curG];
    % comb
    case 'cm'
        % initial gradient locations, will ideally dwindle
        curp = pG;
        % initial gradient values
        curG = G;
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
                curG = [curG(outidx,:); curG(curidx,:)];
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
        gradout = [curp, curG];
end

end