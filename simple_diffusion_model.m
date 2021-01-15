function [z_mod,rss] = simple_diffusion_model(x,z1,z2,k,run_time,varargin)
% simple_diffusion_model will simulate diffusion of a profile given an
% initial and final topographic profile based on an input diffusion 
% parameter, k, and a run_time. The code will evolve the initial
% topographic profile using explicit finite difference calculations and
% compare the results to the final observed profile by calculating the
% residual sum of squares (rss)
%
% Required inputs
% x - distance (in meters) along the profile (not distance should have even increments)
% z1 - elevation (in meters) of the initial topography
% z2 - elevation (in meters) of the final topography
% k - diffusion coefficient in m^2/kyr
% run_time - total model run time in kyrs
%
% Optional inputs
% U - uplift rate in m/kyr; zero is default
% stabil - parameter for numerical stability. if model blows up increase
% value. one is default
% plot_opt - plot results? 0 if no, 1 if yes. 0 is default
%
% Outputs
% z_mod - elevations of the final time step of the model
% rss - residual sum of squares
%
% Examples:
% 
% example 1:
% x = 0:1:100; 
% z1 = ones(1,length(x));
% z1(25:50) = sind(20).*x(1:26)+1;
% z1(51:75) = max(z1) + sind(-20).*x(1:25);
% z2 = ones(1,length(x));
% z2(25:50) = sind(10).*x(1:26)+1;
% z2(51:75) = max(z2) + sind(-10).*x(1:25);
% k = 0.1;
% run_time = 1000;
% [z_mod,rss] = simple_diffusion_model(x,z1,z2,k,run_time);
%
% example 2
% example 1:
% x = 0:1:100; 
% z1 = ones(1,length(x));
% z1(25:50) = sind(20).*x(1:26)+1;
% z1(51:75) = max(z1) + sind(-20).*x(1:25);
% z2 = ones(1,length(x));
% z2(25:50) = sind(10).*x(1:26)+1;
% z2(51:75) = max(z2) + sind(-10).*x(1:25);
% k = 0.1;
% run_time = 1000;
% [z_mod,rss] = simple_diffusion_model(x,z1,z2,k,run_time,'U',0.5,'plot_opt',1);
%
% Author: Sean F. Gallen
% Contact: sean.gallen[at]colostate.edu
% date modified: 10/06/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse the inputs
p = inputParser;
p.FunctionName = 'simple_diffusion_model';

% required inputs
addRequired(p,'x', @(x) isvector(x));
addRequired(p,'z1', @(x) isvector(x));
addRequired(p,'z2', @(x) isvector(x));
addRequired(p,'k', @(x) isscalar(x));
addRequired(p,'run_time', @(x) isscalar(x));

% add optional inputs
addParameter(p,'U', 0, @(x) isscalar(x));
addParameter(p,'stabil', 1, @(x) isscalar(x));
addParameter(p,'plot_opt', 0, @(x) isscalar(x));


% declare variable names
parse(p, x, z1, z2, k, run_time, varargin{:});
x   = p.Results.x;
z1   = p.Results.z1;
z2  = p.Results.z2;
k   = p.Results.k;
run_time = p.Results.run_time;
U   = p.Results.U;
stabil   = p.Results.stabil;
plot_opt  = p.Results.plot_opt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define dx --> assumes dx is uniform
dx = abs(x(2)-x(1));

% model time-step and run length
lambda =0.4/stabil;    % stability parameter that will determine dt based on k and dx.
dt = lambda*dx^2/k;     % calculate dt (n kyrs) to satify the CLF condition
t = [0:dt:run_time]';

% pad the array to avoid edge effects
z_pad = zeros(length(z1)+2,1);
z_pad(1) = z1(1);
z_pad(end) = z1(end);

% set boundary and initial conditions conditions
z_s = z1(1);
z_e = z1(end);
z = z1;

% if we are plotting the data set up the figure
if plot_opt == 1
    figure
    nplots = 10;
    plott = round(length(t)/nplots);
    cols = hsv(nplots);
    legText = cell(1,nplots+2);
    plot(x,z1,'k--','color',[0.5 0.5 0.5],'linewidth',2); hold on
    legText{1} = 'initial profile kyrs';
    n = 1;
end

% run the time look and evolve the profile
for i = 1:length(t)
    
    C = (z_pad(3:end) - 2*z_pad(2:end-1) + z_pad(1:end-2)).*lambda; % calculate curvature of the profile abd diffuse
    z = z + C + U*dt;           % evolve the profile
    z(1) = z_s; z(end) = z_e;   % anchor to boundary conditions
    z_pad(2:end-1) = z;         % reset the padded array
    
    
    % plot results at given time step if needed
    if plot_opt == 1
        if i >= plott*n
            plot(x,z,'-','color',cols(n,:)); hold on
            legText{n+1} = ['time = ', num2str(plott*n*dt),' kyrs'];
            n = n +1;
            pause(0.2)
        end
    end
end

% axes labels and legent if needed
if plot_opt == 1
    plot(x,z2,'k--','linewidth',2);
    if n == 11
        legText{end} = 'final observed profile';
    else
        legText{end-1} = 'final observed profile';
        legText = legText(1:end-1);
    end
    xlabel('distance (m)');
    ylabel('elevation (m)');
    legend(legText);
end

% make the final calculations
z_mod = z;
rss = sum((z2-z_mod).^2);
end