function Z = stratpot2D(grad, trac, P)
%%
% IN: 
% grad: m x 4 array with gradient positions and compoments in columns
% trac: nk x 3 array with trace locations and ids in the columns
% P: l x 2 array of locations for which to interpolate stratigraphic
%    potential
%
% OUT:
% Z: potential interpolated at points P

%% PREPROCESS
% separate positions and measurements
pG = grad(:,1:2);
pZ = trac(:,1:2);

G = grad(:,3:4);
bedID = trac(:,3);

m = size(G,1);
n = nincrements(bedID);

% covariance model, always cubic
range = 1.0;
sill = 0.01;

% drift 2 = linear, 5 = quadratic
drift = 'linear';

% plot or not
plt = 0;

% number of drift functions
switch drift
    case 'linear'
        d = 2;
    case 'quadratic'
        d = 5;
end

%% TEST

% hx = linspace(-2, 2);
% hy = linspace(-2, 2);
% 
% [hX, hY] = meshgrid(hx, hy);
% 
% hx = reshape(hX, numel(hX), 1);
% hy = reshape(hY, numel(hY), 1);
% 
% h = [hx, hy];
% 
% [h_th, h_r] = cart2pol(hx, hy);
% 
% h_pol = [h_th, h_r];
% 
% % plot Cz
% figure
% plot(h_r, Cz(range, sill, h_r))
% hold on
% plot(h_r, dCz(range, sill, h_r))
% plot(h_r, d2Cz(range, sill, h_r))
% legend({'Cz', 'dCz', 'd2Cz'})
% 
% % plot C_GxGy
% CGxGy = C_GxGy(range, sill, h); 
% figure
% contour(hX, hY, reshape(CGxGy, 100, 100))
% title('C_{G_xG_y}')
% colorbar
% 
% % plot C_GxGx
% CGxGx = C_GxGx(range, sill, h); 
% figure
% contour(hX, hY, reshape(CGxGx, 100, 100))
% title('C_{G_xG_x}')
% colorbar

%% COVARIANCES

% gradients
c_gxgx = zeros(m, m);
c_gygy = zeros(m, m);
c_gxgy = zeros(m, m);

c_gxgx_pol = zeros(m, m);
c_gygy_pol = zeros(m, m);
c_gxgy_pol = zeros(m, m);

for ii = 1:m
    for jj = 1:m
        h = pG(ii,:)-pG(jj,:);
        c_gxgx(ii,jj) = C_GxGx(h);
        c_gygy(ii,jj) = C_GyGy(h);
        c_gxgy(ii,jj) = C_GxGy(h);
        
        % do polar coords just for redundancy
        [h_th, h_r] = cart2pol(h(:,1), h(:,2));
        h_pol = [h_th, h_r];
        
        c_gxgx_pol(ii,jj) = C_GxGx_pol(h_pol);
        c_gygy_pol(ii,jj) = C_GyGy_pol(h_pol);
        c_gxgy_pol(ii,jj) = C_GxGy_pol(h_pol);
    end
end

% gradient-increment
c_igx = zeros(n,m);
c_igy = zeros(n,m);

c_igx_pol = zeros(n,m);
c_igy_pol = zeros(n,m);

for ii = 1:n
    
    [rsecondIdx, rfirstIdx] = incrementIdx(bedID, 'increment', ii);
    % cross-covariance, gradient and increments 
    for jj = 1:m
        % distance between gradient and jth point on bed trace
        % need to add directionality here
        h2 = pZ(rsecondIdx,:)-pG(jj,:);
        h1 = pZ(rfirstIdx,:)-pG(jj,:);
        c_igx(ii,jj) = C_IGx(h2) - C_IGx(h1);
        c_igy(ii,jj) = C_IGy(h2) - C_IGy(h1);
        
        % do polar coords just for redundancy
        [h2_th, h2_r] = cart2pol(h2(:,1), h2(:,2));
        h2_pol = [h2_th, h2_r];
        [h1_th, h1_r] = cart2pol(h1(:,1), h1(:,2));
        h1_pol = [h1_th, h1_r];
        
        c_igx_pol(ii,jj) = C_IGx_pol(h2_pol) - ...
                           C_IGx_pol(h1_pol);
        c_igy_pol(ii,jj) = C_IGy_pol(h2_pol) - ...
                           C_IGy_pol(h1_pol);
    end  
   
end

c_ig = [c_igx, c_igy];

% increments
c_i = zeros(n,n);

for ii = 1:n
    
    [rsecondIdx, rfirstIdx] = incrementIdx(bedID, 'increment', ii);
    % cross-covariance, gradient and increments 
    for jj = 1:n
        [ssecondIdx, sfirstIdx] = incrementIdx(bedID, 'increment', jj);
        
        % need total distance, not separation vector, so take norm
        h_rusv = norm(pZ(ssecondIdx,:)-pZ(rsecondIdx,:));
        h_s1ru = norm(pZ(rsecondIdx,:)-pZ(sfirstIdx,:));
        h_svr1 = norm(pZ(rfirstIdx,:)-pZ(ssecondIdx,:));
        h_r1s1 = norm(pZ(sfirstIdx,:)-pZ(rfirstIdx,:));
        
        c_i(ii,jj) = Cz(h_rusv) - ...
                     Cz(h_s1ru) - ...
                     Cz(h_svr1) + ...
                     Cz(h_r1s1);
    end  
   
end

%% DRIFT
% for 2D case, linear drift means basis functions f1 = x, f2 = y
% for 2D case, quadratic drift means basis functions
%   f1 = x, f2 = y, f3 = x^2, f4 = y^2, f5 = x*y
% the constant is already accounted for by assuming an origin value of 0

% gradients
switch drift
    case 'linear'
        u_gx = [ones(1,m); zeros(1,m)];
        u_gy = [zeros(1,m); ones(1,m)];
        
    case 'quadratic'
        u_gx = [ones(1,m); 
                zeros(1,m);
                2*pG(:,1)';
                zeros(1,m);
                pG(:,2)'];
        u_gy = [zeros(1,m); 
                ones(1,m);
                zeros(1,m);
                2*pG(:,2)';
                pG(:,1)'];
end

u_g = [u_gx, u_gy];


% interfaces
switch drift
    case 'linear'
        u_i = zeros(2, n);
    case 'quadratic'
        u_i = zeros(5, n);
end

for ii = 1:n
    [secondIdx, firstIdx] = incrementIdx(bedID, 'increment', ii);
    
    switch drift
        case 'linear'
            u_i(1,ii) = pZ(secondIdx,1) - pZ(firstIdx,1);
            u_i(2,ii) = pZ(secondIdx,2) - pZ(firstIdx,2);
            
        case 'quadratic'
            u_i(1,ii) = pZ(secondIdx,1) - pZ(firstIdx,1);
            u_i(2,ii) = pZ(secondIdx,2) - pZ(firstIdx,2);
            u_i(3,ii) = pZ(secondIdx,1)^2 - pZ(firstIdx,1)^2;
            u_i(4,ii) = pZ(secondIdx,2)^2 - pZ(firstIdx,2)^2;
            u_i(5,ii) = pZ(secondIdx,1)*pZ(secondIdx,2) - ...
                        pZ(firstIdx,1)*pZ(firstIdx,2);
            
    end
end

u = [u_g, u_i];

%% SET UP AND SOLVE COKRIGING SYSTEM

% covariances
c_g = [c_gxgx, c_gxgy';
       c_gxgy, c_gygy];
c = [c_g, c_ig';
     c_ig, c_i];

% cokriging matrix
K = [c, u';
     u, zeros(d,d)];

% right hand side, "data"
D = [G(:,1); G(:,2); zeros(n+d, 1)];

% get weights
w = K\D;

% split 
w_Gx = w(1:m);
w_Gy = w(m+1:2*m);
w_i = w(2*m+1:2*m+n);
w_d = w(2*m+n+1:end);


%% EVALUATE ON GRID

ninterp = size(P,1);
Z = zeros(ninterp, 1);

for ii = 1:ninterp
    cur_p = P(ii, :);
    % gradient contributions
    for jj = 1:m
        cur_pG = pG(jj, :);
        Z(ii) = Z(ii) + w_Gx(jj)*C_IGx(cur_p-cur_pG) + ...
                        w_Gy(jj)*C_IGy(cur_p-cur_pG);
    end
    % increment contributions
    for jj = 1:n
        [secondIdx, firstIdx] = incrementIdx(bedID, 'increment', jj);
        Z(ii) = Z(ii) + w_i(jj)*(Cz(norm(cur_p-pZ(secondIdx,:))) - ...
                                 Cz(norm(cur_p-pZ(firstIdx,:))));
    end
    % drift contributions
    switch drift
        case 'linear'
            Z(ii) = Z(ii) + w_d(1)*cur_p(1) + w_d(2)*cur_p(2);
        case 'quadratic'
            Z(ii) = Z(ii) + w_d(1)*cur_p(1) + w_d(2)*cur_p(2) + ...
                            w_d(3)*cur_p(1)^2 + w_d(4)*cur_p(2)^2 + ...
                            w_d(5)*cur_p(1)*cur_p(2);
    end
end

%% VISUALIZE

if plt
    
    X = reshape(X,nint,nint);
    Y = reshape(Y,nint,nint);
    Zgrd = reshape(Z,nint,nint);
    [dxZgrd,dyZgrd] = gradient(Zgrd);

    col = linspecer(length(unique(bedID)));

    figure
    contourf(X,Y,Zgrd,30)
    hold on
    quiver(pG(:,1),pG(:,2),G(:,1),G(:,2),0.5,'k')
    quiver(X,Y,dxZgrd,dyZgrd,2)
    scatter(pZ(:,1),pZ(:,2), 50, col(bedID,:), 'filled','MarkerEdgeColor', 'k')
    xlabel('x')
    ylabel('y')
    axis equal
    
end

%% POTENTIAL COVARIANCE FUNCTION AND DERIVATIVES

function C = Cz(r)
    
    idx = r >= range;
    C(~idx) = sill * ...
                (1 - ...
                7*(r(~idx)/range).^2 + ...
                (35/4)*(r(~idx)/range).^3 - ...
                (7/2)*(r(~idx)/range).^5 + ...
                (3/4)*(r(~idx)/range).^7 );
    C(idx) = 0;
    
    % ensure column
    C = C(:);
    
end

% first derivative of covariance with respect to r, multiplied by 1/r since
% it always appears as 1/r * dCz
function C = dCz(r)
    idx = r >= range;
    C(~idx) = sill/(range^2) .* ...
        (-14 + (105/4)*(r(~idx)/range) - ...
                (35/2)*(r(~idx)/range).^3 + ...
                (21/4)*(r(~idx)/range).^5 );
    C(idx) = 0;
    
    % ensure column
    C = C(:);
end

% second derivative of covariance with respect to r
function C = d2Cz(r)

    idx = r >= range;
    C(~idx) = sill/(range^2) * ...
        (-14 + (105/2)*(r(~idx)/range) - ...
                    70*(r(~idx)/range).^3 + ...
                (63/2)*(r(~idx)/range).^5 );
    C(idx) = 0;
    
    % ensure column
    C = C(:);

end


%% GRADIENT COVARIANCE AND CROSS COVARIANCE FUNCTIONS
% recall that dCz is already multiplied by 1/r


% x-component covariance
% h is nx2 array with first column theta and second column r
function C = C_GxGx_pol(h)
    
    C = cos(h(:,1)).^2 .* (dCz(h(:,2)) - ...
                           d2Cz(h(:,2))) - ...
        dCz(h(:,2));
        
end

% y-component covariance
% h is nx2 array with first column theta and second column r
function C = C_GyGy_pol(h)
    
    C = sin(h(:,1)).^2 .* (dCz(h(:,2)) - ...
                           d2Cz(h(:,2))) - ...
        dCz(h(:,2));
        
end

% gradient x-y cross covariance
% h is nx2 array with first column theta and second column r
function C = C_GxGy_pol(h)
    
    C = cos(h(:,1)) .* sin(h(:,1)) .* ...
        ( dCz(h(:,2)) - ...
          d2Cz(h(:,2)) ); 
        
end

% x-component covariance
% h is nx2 array with first column hx and second column hy
function C = C_GxGx(h)
    
    nh = size(h, 1);
    C = zeros(nh, 1);
    
    r = sqrt(sum(h.^2, 2));
    idx = r ~= 0;
    
    C(~idx) = 14*sill/(range^2);
    
    C(idx) = (h(idx,1).^2./r(idx).^2) .* (dCz(r(idx)) - d2Cz(r(idx))) - ...
             dCz(r(idx));
        
end

% y-component covariance
% h is nx2 array with first column hx and second column hy
function C = C_GyGy(h)
    
    nh = size(h, 1);
    C = zeros(nh, 1);
    
    r = sqrt(sum(h.^2, 2));
    idx = r ~= 0;
    
    C(~idx) = 14*sill/(range^2);
    
    C(idx) = (h(idx,2).^2./r(idx).^2) .* ...
                (dCz(r(idx)) - d2Cz(r(idx))) - ...
             dCz(r(idx));
       
end

% gradient x-y cross covariance
% h is nx2 array with first column hx and second column hy
function C = C_GxGy(h)
    
    nh = size(h, 1);
    C = zeros(nh, 1);
    
    r = sqrt(sum(h.^2, 2));
    idx = r ~= 0;
    
    C(~idx) = 0;
    
    C(idx) = (h(idx,1).*h(idx,2)./r(idx).^2) .* ...
                (dCz(r(idx)) - d2Cz(r(idx)));
 
end

%% GRADIENT-INCREMENT COVARIANCE
% recall that dCz is already multiplied by 1/r

% increment x-gradient covariance
% h is nx2 array with first column theta and second column r
function C = C_IGx_pol(h)

    C = -h(:,2) .* cos(h(:,1)) .* dCz(h(:,2));
    
end

% increment y-gradient covariance
% h is nx2 array with first column theta and second column r
function C = C_IGy_pol(h)

    C = -h(:,2) .* sin(h(:,1)) .* dCz(h(:,2));
    
end

% increment x-gradient covariance
% h is nx2 array with first column hx and second column hy
function C = C_IGx(h)

    r = sqrt(sum(h.^2, 2));
    C = -h(:,1) .* dCz(r);
    
end

% increment y-gradient covariance
% h is nx2 array with first column hx and second column hy
function C = C_IGy(h)

    r = sqrt(sum(h.^2, 2));
    C = -h(:,2) .* dCz(r);
    
end

end