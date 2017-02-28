%% DUAL FORM 
% example is in dual form
function stratpotexample2D()


grad = load('gradient.csv');
trac = load('trace.csv');

%% PREPROCESS
% separate positions and measurements
pG = grad(:,1:2);
pZ = trac(:,1:2);

G = grad(:,3:4);
bedID = trac(:,3);

m = size(G,1);
n = nincrements(bedID);

% covariance model
range = 1;
sill = 1;
mod = 'cubic';  % Cubic or Gaussian

% drift 2 = linear, 5 = quadratic
drf = 5;

%% DRIFT BASIS FUNCTIONS

if drf == 2
    % linear
    f{1} = @(x,y) x;
    f{2} = @(x,y) y;
    dfx{1} = @(x,y) 1;
    dfx{2} = @(x,y) 0;
    dfy{1} = @(x,y) 0;
    dfy{2} = @(x,y) 1;
    ord = 'linear';
elseif drf == 5
    % quadratic
    f{1} = @(x,y) x;
    f{2} = @(x,y) y;
    f{3} = @(x,y) x.^2;
    f{4} = @(x,y) y.^2;
    f{5} = @(x,y) x.*y;
    dfx{1} = @(x,y) 1;
    dfx{2} = @(x,y) 0;
    dfx{3} = @(x,y) 2*x;
    dfx{4} = @(x,y) 0;
    dfx{5} = @(x,y) y;
    dfy{1} = @(x,y) 0;
    dfy{2} = @(x,y) 1;
    dfy{3} = @(x,y) 0;
    dfy{4} = @(x,y) 2*y;
    dfy{5} = @(x,y) x;
    ord = 'quadratic';
end

%% GENERATE KRIGING MATRIX

% construct covariance matrices from functions
Cgxx = zeros(m,m);
Cgyy = zeros(m,m);
Cgxy = zeros(m,m);
for j = 1:m
    % gradients
    for k = 1:m
        dist = (pG(j,:)-pG(k,:));
        Cgxx(j,k) = cgxx(range,sill,dist);
        Cgyy(j,k) = cgyy(range,sill,dist);
        Cgxy(j,k) = cgxy(range,sill,dist);
    end
end
Czgx = zeros(n,m);
Czgy = zeros(n,m);
Cz = zeros(n,n);
for j = 1:n
    [rsecondIdx, rfirstIdx] = incrementIdx(bedID,j);
    % cross-covariance, gradient and increments 
    for k = 1:m
        % distance between gradient and jth point on bed trace
        % need to add directionality here
        dist2 = pZ(rsecondIdx,:)-pG(k,:);
        dist1 = pZ(rfirstIdx,:)-pG(k,:);
        Czgx(j,k) = czgx(range,sill,dist2) - czgx(range,sill,dist1);
        Czgy(j,k) = czgy(range,sill,dist2) - czgy(range,sill,dist1);
    end  
    % increments
    for k = 1:n
        [ssecondIdx, sfirstIdx] = incrementIdx(bedID,k);
        % take norms for cz, which wants only r
        distrsjk = norm(pZ(rsecondIdx,:)-pZ(ssecondIdx,:));
        distrs1k = norm(pZ(sfirstIdx,:)-pZ(rsecondIdx,:));
        distrsj1 = norm(pZ(ssecondIdx,:)-pZ(rfirstIdx,:));
        distrs11 = norm(pZ(rfirstIdx,:)-pZ(sfirstIdx,:));
        Cz(j,k) = cz(range,sill,distrsjk) - cz(range,sill,distrs1k) - ...
            cz(range,sill,distrsj1) + cz(range,sill,distrs11);
    end
    
end


% drift function
Fg = gradientDrift(pG,'order',ord);
Fi = incrementDrift(pZ,bedID,'order',ord);
F = [Fg, Fi];

% final kriging matrix
CG = [Cgxx Cgxy'; Cgxy, Cgyy];
CGI = [Czgx, Czgy];
C = [CG, CGI'; CGI, Cz];
K = [C, F'; F, zeros(drf,drf)];

%% SOLVE SYSTEM FOR WEIGHTS

B = [G(:,1); G(:,2); zeros(n,1); zeros(drf,1)];
w = K\B;

a = w(1:m);         % weights on Czgx
b = w(m+1:2*m);     % weights on Czgy
c = w(2*m+1:end-drf); % weights on increments
d = w(end-drf+1:end); % weights on basis functions

%% GRID A SOLUTION SPACE TO INTERPOLATE AT

clear x y X Y
nint = 50;
x = linspace(-0.2,1.2,nint);
y = x;

[X,Y] = meshgrid(x,y);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);

P = [X,Y];

%% INTERPOLATE AT REQUESTED POINTS

% interpolate potential on grid
Zgrd = krigeZ(P);
pZval = krigeZ(pZ);

% interpolate gradient
Gxval = krigeDZx(pG);
Gyval = krigeDZy(pG);

%% RESHAPE AND PLOT

X = reshape(X,nint,nint);
Y = reshape(Y,nint,nint);
Zgrd = reshape(Zgrd,nint,nint);
[dxZgrd,dyZgrd] = gradient(Zgrd);

figure
contourf(X,Y,Zgrd,30)
hold on
quiver(pG(:,1),pG(:,2),G(:,1),G(:,2),0.5,'k')
quiver(X,Y,dxZgrd,dyZgrd,2)
plot(pZ(:,1),pZ(:,2),'ko','markerfacecolor','k')
xlabel('x')
ylabel('y')
axis equal

%% AUXILIARY FUNCTIONS

%% kriging functions
    % krige potential
    function  Z = krigeZ(p)
        Z = zeros(size(p,1),1);
        for jj = 1:size(p,1)
            % four terms in interpolator
            terms = zeros(4,1);
            % gradient-increment cross-covariance
            for kk = 1:m
                distfun = (p(jj,:) - pG(kk,:));
                terms(1) = terms(1) + ...
                    a(kk)*czgx(range,sill,distfun);
                terms(2) = terms(2) + ...
                    b(kk)*czgy(range,sill,distfun) ;
            end
            % increment covariance
            for kk = 1:n
                [secondIdx,firstIdx] = incrementIdx(bedID,kk);
                % must take norm for cz
                dist2fun = norm(p(jj,:)-pZ(secondIdx,:));
                dist1fun = norm(p(jj,:)-pZ(firstIdx,:));
                terms(3) = terms(3) + ...
                   c(kk)* ( cz(range,sill,dist2fun) - ...
                              cz(range,sill,dist1fun) );
            end
            % drift
            for kk = 1:drf
                terms(4) = terms(4) + d(kk)*f{kk}(p(jj,1),p(jj,2));
            end
            Z(jj) = sum(terms);
        end
    end
    % krige x component of gradient of potential
    function dZx = krigeDZx(p)
        dZx = zeros(size(p,1),1);
        for jj = 1:size(p,1)
            % four terms in interpolator
            terms = zeros(4,1);
            % gradient-increment cross-covariance
            for kk = 1:m
                distf = p(jj,:) - pG(kk,:);
                terms(1) = terms(1) + ...
                    a(kk)*cgxx(range,sill,distf);
                terms(2) = terms(2) + ...
                    b(kk)*cgxy(range,sill,distf);  
            end
            % increment covariance
            for kk = 1:n
                [secondIdx,firstIdx] = incrementIdx(bedID,kk);
                dist2f = p(jj,:)-pZ(secondIdx,:);
                dist1f = p(jj,:)-pZ(firstIdx,:);
                terms(3) = terms(3) + ...
                    c(kk)*( czgx(range,sill,dist2f) - czgx(range,sill,dist1f) );
            end
            % drift
            for kk = 1:drf
                terms(4) = terms(4) + d(kk)*dfx{kk}(p(jj,1),p(jj,2));
            end
            dZx(jj,1) = sum(terms);
        end
    end
    % krige y component of gradient of potential
    function dZy = krigeDZy(p)
        dZy = zeros(size(p,1),1);
        for jj = 1:size(p,1)
            % four terms in interpolator
            terms = zeros(4,1);
            % gradient-increment cross-covarianc
            for kk = 1:m
                distf = p(jj,:) - pG(kk,:);
                terms(1) = terms(1) + ...
                    a(kk)*cgxy(range,sill,distf);
                terms(2) = terms(2) + ...
                    b(kk)*cgyy(range,sill,distf) ;
            end
            % increment covariance
            for kk = 1:n
                [secondIdx,firstIdx] = incrementIdx(bedID,kk);
                dist2f = p(jj,:)-pZ(secondIdx,:);
                dist1f = p(jj,:)-pZ(firstIdx,:);
                terms(3) = terms(3) + ...
                    c(kk)*( czgy(range,sill,dist2f) - czgy(range,sill,dist1f) );
            end
            % drift
            for kk = 1:drf
                terms(4) = terms(4) + d(kk)*dfy{kk}(p(jj,1),p(jj,2));
            end
            dZy(jj,1) = sum(terms);
        end
    end

%% covariance functions

% base covariance function (potential covariance, isotropic)
    function C = cz(range,sill,r)
        switch mod
            case 'cubic'
                idx = r >= range;
                C(~idx) = sill * (1 - 7*(r/range)^2 + ...
                           (35/4)*(r/range)^3 - ...
                            (7/2)*(r/range)^5 + ...
                            (3/4)*(r/range)^7);
                C(idx) = 0;
            case 'gaussian'
                C = sill*exp(-norm(r)^2/range^2);
            otherwise
                error('chosen model does not exist')
        end
    end

% first derivative of covariance with respect to r
    function C = dcz(range,sill,r)
        switch mod
            case 'cubic'
                idx = r >= range;
                C(~idx) = sill/range^2 * r * ...
                    (-14 + (105/4)*(r/range) - ...
                            (35/2)*(r/range)^3 + ...
                            (21/4)*(r/range)^5 );
                C(idx) = 0;
            case 'gaussian'
        
        end
    end

% second derivative of covariance with respect to r
    function C = d2cz(range,sill,r)
        switch mod
            case 'cubic'
                idx = r >= range;
                C(~idx) = sill/range^2 * ...
                    (-14 + (105/2)*(r/range) - ...
                                70*(r/range)^3 + ...
                            (63/2)*(r/range)^5 );
                C(idx) = 0;
            case 'gaussian'
        
        end
    end

% gradient covariance functions
    function C = cgxx(range,sill,h)
        switch mod
            case 'cubic'
                if norm(h) == 0
                    C = 14*sill/range^2;
                else
                    r = norm(h);
                    C = (h(1)^2/r^3 - 1/r)*dcz(range,sill,r) - ...
                        (h(1)/r)^2 * d2cz(range,sill,r);
                end
            case 'gaussian'
                C = ( 2/range^2 - 4*h(1)^2/range^4 ) * cz(range,sill,h);
            otherwise
                error('chosen model does not exist')
        end
    end
    function C = cgyy(range,sill,h)
        switch mod
            case 'cubic'
                if norm(h) == 0
                    C = 14*sill/range^2;
                else
                    r = norm(h);
                    C = (h(2)^2/r^3 - 1/r)*dcz(range,sill,r) - ...
                        (h(2)/r)^2 * d2cz(range,sill,r);
                end
            case 'gaussian'
                C = ( 2/range^2 - 4*h(2)^2/range^4 ) * cz(range,sill,h);
            otherwise
                error('chosen model does not exist')
        end
    end
    function C = cgxy(range,sill,h)
        switch mod
            case 'cubic'
                % will need to change this for 3-d case?
                if norm(h) == 0
                    C = 0;
                else
                    r = norm(h);
                    C = h(1)*h(2)/r^2 * ...
                        ((1/r)*dcz(range,sill,r) - ...
                            d2cz(range,sill,r));
                end
            case 'gaussian'
                C = ( -4*h(1)*h(2)/range^4 ) * cz(range,sill,h);
            otherwise
                error('chosen model does not exist')
        end
    end
% gradient-potential cross-covariance models
    function C = czgx(range,sill,h)
        switch mod
            case 'cubic'
                r = norm(h);
                C = -h(1)/r * dcz(range,sill,r);
            case 'gaussian'
                C = ( 2*h(1)/range^2 ) * cz(range,sill,h);
            otherwise
                error('chosen model does not exist')
        end
    end
    function C = czgy(range,sill,h)
        switch mod
            case 'cubic'
                r = norm(h);
                C = -h(2)/r * dcz(range,sill,r);
            case 'gaussian'
                C = ( 2*h(2)/range^2 ) * cz(range,sill,h);
            otherwise
                error('chosen model does not exist')
        end
    end

end