clear
clc

% load the excel file
data = readtable('profile_data.xlsx');

% unpack the data
x = data.distance_m_;           % distance along the profile
zi = data.initialElevation_m_;  % initial elevation of the observed profile
zf = data.finalElevation_m_;    % final elevation of the observed profile

% define parameters
k = 0.1;            % diffusion parameter m^2/kyr
run_time = 50;      % model run time in kyr

% run the model
[z_mod,rss] = simple_diffusion_model(x,zi,zf,k,run_time,'plot_opt',1);



