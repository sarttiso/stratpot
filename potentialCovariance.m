% this function is responsible for producing functional forms of nested
% potential covariance models based on given ranges and sills. The ranges
% and sills reflect ranges and sills in given directions (when considering
% anisotropy) for each gradient component; these can be found with the
% variogram_fit function. 
%
% range_all: ndir x 3 matrix where each row corresponds to ranges for each
%   gradient component [Gx,Gy,Gz] in a given direction.
% sill_all: ndir x 3 matrix where each row corresponds to sills for each
%   gradient component [Gx,Gy,Gz] in a given direction.
% 'derivative': parameter specifying which derivative of the potential
%   covariance function to return (default 0 is the base covariance
%   function), only first and second derivatives are supported
% 'model': what potential covariance model to use (default 'cubic')
%
% Cnest: matrix of nested potential covariance functions, each
%   corresponding to the directions corresponding to the input rows of
%   range_all and sill_all.
% 
% LIMITATIONS: 
%   - currently only supports the cubic model as in Aug (2004,2005)

function Cnest = potentialCovariance(range_all,sill_all,varargin)
    parser = inputParser;
    
    addRequired(parser,'range_all',@isnumeric)
    addRequired(parser,'sill_all',@isnumeric)
    addParameter(parser,'derivative',0,@(x) any(x == [0,1,2]))
    addParameter(parser,'model','cubic',@ischar)
%     addParameter(parser,'ndim',3,@isnumeric)
    
    parse(parser,range_all,sill_all,varargin{:});
    
    range_all = parser.Results.range_all;
    sill_all = parser.Results.sill_all;
    derivative = parser.Results.derivative;
    model = parser.Results.model;
%     ndim = parser.Results.ndim;
    
    assert(all(size(range_all) == size(sill_all)), ...
        'ranges and sills are not the same size')
    assert(size(range_all,2) == 3 && size(sill_all,2) == 3, ...
        'ranges and sills must have 3 columns')
    
    model = validatestring(model,{'cubic'});
    ndir = size(range_all,1);
    
    % get base potential covariance model to combine in a nested structure
    if strcmp(model,'cubic')
        switch derivative
            % base cubic function
            case 0
                cz = @cubic;
            % first derivative
            case 1 
                cz = @dcubic;
            % second derivative
            case 2
                cz = @d2cubic;                             
        end
    end
    
    % create nested structure for each direction
    for j = 1:ndir
        % order gradient component sills from smallest to largest
        ord = sortrows([sill_all(j,:); 1 2 3;]',1);
        % get indices of increasing sills
        ordidx = ord(:,2);
        % get ordered ranges
        range_pot = range_all(j,ordidx);
        % convert gradient component sills to potential sills
        %   first get increments in sill for nested components such that the
        %   sill of the largest component is the sum of the sills of the
        %   nested components.
        sillinc = diff(sill_all(j,ordidx)); 
        sillinc = [sill_all(j,ordidx(1)) sillinc];
        %   then convert to potential sills using corresponding ranges,
        %   which depends also on the selected covariance model
        if strcmp(model,'cubic')
            sill_pot = (sillinc.*range_pot.^2)/14;
        end
        
        % h is a 3-vector
        Cnest{j} = @(h) ...
         cz(range_pot(1),sill_pot(1),sqrt(sum(h.^2,2))) + ...          % K3
         cz(range_pot(2),sill_pot(2),sqrt(sum(h(:,ordidx(2:3)).^2,2)))+...% K2
         cz(range_pot(3),sill_pot(3),abs(h(:,ordidx(3))));            % K1
    end

%%  AUXILIARY FUNCTIONS
% covariance functions
    % base covariance function for cubic case
    function cz = cubic(a,c0,r)
        idx = r >= a;
        cz(idx) = 0;
        cz(~idx) = c0 * ...
            (1 - 7*(r(~idx)/a).^2 + (35/4)*(r(~idx)/a).^3 - ...
             (7/2)*(r(~idx)/a).^5 + (3/4)*(r(~idx)/a).^7);
        cz = cz';
    end
    % first derivative of the potential covariance function, BUT normalized
    % by 1/r, since all of the derived cross-covariance terms depend always
    % on 1/r dC/dr
    function cz = dcubic(a,c0,r)
        idx = r >= a;
        cz(idx) = 0;
        cz(~idx) = c0/a^2 * ...
            (-14 + (105/4)*(r(~idx)/a) - ...
            (35/2)*(r(~idx)/a).^3 + ...
            (21/4)*(r(~idx)/a).^5 );
        cz = cz';
    end
    % second derivative of cubic covariance function, no normalization
    function cz = d2cubic(a,c0,r)
        idx = r >= a;
        cz(idx) = 0;
        cz(~idx) = c0/a^2 * ...
            (-14 + (105/2)*(r(~idx)/a) - ...
                        70*(r(~idx)/a).^3 + ...
                    (63/2)*(r(~idx)/a).^5 );
        cz = cz';
    end
    
    
end