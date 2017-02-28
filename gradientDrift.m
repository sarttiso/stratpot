% given a set of positions in an input vector p, generate the matrix
% corresponding to the evaluation of the drift basis functions at gradient
% data locations.
% p: npositions (rows) x ncoordinates (columns)
% 'order' : 'linear' or 'quadratic', defaults to 'linear'
% LIMITATIONS:
%   - currently written for dual form; might need modification for primal
%   form (i.e. column of ones corresponding to constant drift, e.g. mean)

function Fg = gradientDrift(p,varargin)
    parser = inputParser;
    addRequired(parser,'p',@isnumeric);
    addParameter(parser,'order','linear',@ischar);
    
    parse(parser,p,varargin{:});
    
    % gradient data locations 
    p = parser.Results.p;
    m = size(p,1);  % number of input data
    D = size(p,2);  % dimension of input data space
    
    % drift function order 
    order = parser.Results.order;
    order = validatestring(order,{'linear','quadratic'});
    
    if strcmp(order,'linear')
        Fg = zeros(D,D*m);
        for j = 1:D
            Fg(j,m*(j-1)+1:m*j) = ones(1,m);
        end
        
    elseif strcmp(order,'quadratic') && D == 2
        Fg = zeros(5,D*m);
        % fill rows 1-4
        for j = 1:2
            Fg(j,m*(j-1)+1:m*j) = ones(1,m);
            Fg(j+2,m*(j-1)+1:m*j) = 2*p(:,j);
        end
        % fill 5th row
        Fg(5,1:m) = p(:,2);
        Fg(5,m+1:end) = p(:,1);
        
    elseif strcmp(order,'quadratic') && D == 3
        Fg = zeros(9,D*m);
        % fill rows 1-6
        for j = 1:3
            Fg(j,m*(j-1)+1:m*j) = ones(1,m);
            Fg(j+3,m*(j-1)+1:m*j) = 2*p(:,j);
        end
        % fill rows 7-9
        Fg(7,1:m)       = p(:,2);
        Fg(7,m+1:2*m)   = p(:,1);
        Fg(8,1:m)       = p(:,3);
        Fg(8,2*m+1:end) = p(:,1);
        Fg(9,m+1:2*m)   = p(:,3);
        Fg(9,2*m+1:end) = p(:,2);
    else
        error('check input data')
    end
    
end