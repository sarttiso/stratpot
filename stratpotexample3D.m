% function stratpotexample3D()

grad = load('gradientgrid3D.csv');
trac = load('tracegrid3D.csv');
    
%% PREPROCESS
% separate positions and measurements
pG = grad(:,1:3);
pZ = trac(:,1:3);

G = grad(:,4:6);
G = bsxfun(@rdivide,G,sqrt(sum(G.^2,2)));
bedID = trac(:,4);

% covariance model
range = 2;
sill = 1e-6;

% drift linear, quadratic
drf = 'linear';

% points at which to interpolate
n = 100;
x = linspace(0,1,n);
y = x;
z = x;
[X,Y,Z] = meshgrid(x,y,z);
x = reshape(X,numel(X),1);
y = reshape(Y,numel(Y),1);
z = reshape(Z,numel(Z),1);
pS = [x,y,z];

pot = stratpot(pZ,bedID,pG,G,pS,range,sill, ...
    'usevarioparams',true,'krigdrift',drf);

% pot = reshape(pot,n,n,n);
% [dxpot,dypot,dzpot] = gradient(pot);


% end