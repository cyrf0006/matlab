%%  Kelvin Helmholtz instability
%  Two layers with uniform and opposite velocities (sheared flow) are
% separated by a transition zone where the density linearly and the
% velocity vary linearly.  User can change the width of the transition
% zone and the Richardson number
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Author: G. Roullet (roullet@univ-brest.fr)
% Creation date: October 2013
%
clear param

% physics
param.model = 'nh';
param.noslip = 0; % 1=activate no slip BC (boundary layers)
		  % 0=set a free slip BC (inviscid case)

% physical constants
param.g  = 1; % gravity

% domain and resolution      
param.nx = 300;
param.ny = param.nx/2;
param.Lx = 8; % domain size
param.Ly = 4;


% time
param.tend = 60;

% plotting options
param.nplot = 8; % plotting frequency  (one every nplot time steps)
param.varplot = 5; % index of variable to plot 
                   % U=1 / V=2 / vorticity=3 / PSI=4 / buoyancy=5
param.color=2;
param.cax=[-1. 1.]*.1;
param.typeplot='pcolor';
param.plotpsi=0;
param.movie=0;



% numerics
param.cholesky=1; % use a cholesky factorization 
                  % to speed up the vorticity inversion
param.timestepping = 'RK3'; % 'Heun', 'RK2', 'RK3'
param.upwind = 'up5'; % 'up1' or 'up3' or 'up5' 
param.splitting='lxf';% 'lxf' or 'max'
param.fixed_dt=1;
param.cfl = 0.8; % desired CFL, used if param.fixed_dt=0
param.dt=.0125; % used if param.fixed_dt=1

%% display param on the command window
param

%% grid
dx = param.Lx/param.nx;
dy = param.Ly/param.ny;
x = (0.5:param.nx-.5)*dx;
y = (0.5:param.ny-.5)*dy;
[xr,yr] = meshgrid(x,y); % T-point

param.dx = dx;
param.dy = dy;
%% set solid boundaries
global index

% index is defined on T points
% it assigns each solid object an index

index=zeros(param.ny,param.nx);

index(1,:) = 1; % periodic domain
index(end,:) = 2; 


% set circulations around objects
param.circ=[-1 1];



%% deduce masks from index array
% and set up the state vector

init

%% define a BW colormap switching from black to white with a narrow
% gray transition
figure(1)
cmap=colormap;
ncol=size(cmap,1);
ii=1:ncol;
cmap(:,1)= 0.5*(1+tanh( (ii-1-(ncol-1)/2)/(ncol-1)*5 ));
cmap(:,2)=cmap(:,1);
cmap(:,3)=cmap(:,1);
colormap(cmap)

g=gray(64);
gg=flipud(g);
h=hot(64);
colormap([g;gg(1:8:end,:);h])


%% initial state [mask can be used here]
h = zeros(param.ny,param.nx);

% two layer fluid + initial perturbation
amp=1e-2;
%z=yr-param.Ly/2+ param.Ly * amp * cos(8*pi*xr/param.Lx);
%z=yr-param.Ly*(.1+amp/2)+param.Ly/2 * amp * tanh( (xr/param.Lx-.5)/.1);

z=yr-param.Ly/2-param.Ly*amp*sin(xr*2*pi/param.Lx);

% discontinuity thickness
param.thickness=0.2/param.dx;% in number of points
delta_u=1;
delta_b=.1;
richardson = param.g*delta_b*(param.thickness*param.dx)/(delta_u^2);

richardson = .15;
delta_b = richardson*(delta_u^2)/(param.g*param.thickness*param.dx)
%z=yr-param.Ly/2- param.Ly * amp * exp( -(xr/param.Lx-.5).^2 / (2*.04^2) );

h=tanh( (z/param.Ly)*param.ny/param.thickness);


% velocity jump is +2 (from -1 to 1)
u=h*delta_u;
v=h*0;
h=h*delta_b;

% rescale coloraxis for buoyancy
param.cax=[-1. 1.]*delta_b;


param.circ=[1 -1]*param.nx*delta_u;

vor=(GTy'*u(:)-GTx'*v(:));

s(:,param.vorticity) = vor(:);


% half buoyancy jump


%% ok that's it, we have enough information to continue

multimod

%% analysis