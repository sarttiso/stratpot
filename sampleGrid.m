% this function generates a sample grid 

function sampleGrid(xmin,xmax,ymin,ymax,varargin)
    parser = inputParser;
    addRequired(parser,'xmin',@isnumeric);
    addRequired(parser,'xmax',@isnumeric);
    addRequired(parser,'ymin',@isnumeric);
    addRequired(parser,'ymax',@isnumeric);
    addOptional(parser,'dx',1,@isnumeric);
    addOptional(parser,'dy',[],@isnumeric);
    
    parse(parser,xmin,xmax,ymin,ymax,varargin{:});
    
    xmin = parser.Results.xmin;
    xmax = parser.Results.xmax;
    ymin = parser.Results.ymin;
    ymax = parser.Results.ymax;
    dx = parser.Results.dx;
    dy = parser.Results.dy;
    
    assert(xmin < xmax, 'xmin must be less than xmax')
    assert(ymin < ymax, 'xmin must be less than xmax')
    
    if isempty(dy)
        dy = dx;
    end
    assert(dx < xmax-xmin, 'dx cannot be greater than xlimits')
    assert(dy < ymax-ymin, 'dy cannot be greater than ylimits')
end