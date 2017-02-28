function O = gxx_test(b,r,varargin)
    parser = inputParser;
    addRequired(parser,'b',@isnumeric);
    addRequired(parser,'r',@isnumeric);
    addOptional(parser,'theta',0);
    
    parse(parser,b,r,varargin{:});
    
    b = parser.Results.b;
    r = parser.Results.r;
    theta = parser.Results.theta;
    
    O = b(2) * ...
        (cos(theta)^2*(15/4*r./b(1)-5*(r./b(1)).^3+9/4*(r./b(1)).^5) + ...
         sin(theta)^2*(15/8*r./b(1)-5/4*(r./b(1)).^3+3/8*(r./b(1)).^5));
end