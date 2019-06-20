%% Relative Pose Tracker Script
% Author:   Siddarth Kaki

%% housekeeping
clear; close all; clc;

%% options
options.plot = 0;
options.generate_measurements = 1;
options.write_to_file = 1;
options.time = 1;
%#ok<*NOCOL>
   
%% params
% camera params
params.cam.f = 35*10^(-3); % camera focal length

% LM params
params.lm.lambda = 5;
params.lm.max_count = 10000;
params.lm.eps = 0.01;
params.lm.num_init = 10;
params.lm.reinit_att_noise_std = 2;

%% init
rng(2)

%% construct feature points
% specify rigid position vector of feature points wrt target in target frame
rFeaMat = [ 0, 0, 0.5;
            0, 0,-1.5;
            0, 1, 1;
            0,-1, 1; ];

%rFeaMat = [rFeaMat; [1, 0, 0.5]]; % additional feature points
        
[numPts,~] = size(rFeaMat); 

%% setup relative geometry
% specify rigid position vector of camera wrt chaser in chaser frame
rCamVec = [0, 0, 0]';

%% specify true relative state
% rRelVec = [0, 0, 20]'; % Note: z-dir is camera boresight in chaser frame
% phi = pi/6; % roll
% theta = pi/6; % pitch
% psi = pi/6; % yaw
% xVec = [rRelVec; phi; theta; psi];
% xMat = xVec;

filename_read = 'data/pose_trajectory_true.txt';
xMat = read_poses(filename_read,1,linecount(filename_read)); % read in pose trajectory to run

if options.write_to_file == 1,
    % setup output file to write to
    delete data/pose_trajectory_estimated.txt
    fileID_write = fopen('data/pose_trajectory_estimated.txt','a');
    
    delete data/conj_pose_trajectory_estimated.txt
    fileID_write2 = fopen('data/conj_pose_trajectory_estimated.txt','a');
    fmt_write = '%d %d %d %d %d %d\n';
end

if options.time == 1, tic; end

xHatMat = zeros(size(xMat));
xHatMat(1,:) = [0; 0; 30; 0; 0; 0]';

%% loop over poses
for idx = 1:size(xMat,1),
    
    % extract i-th pose
    xVec = xMat(idx,:)';

    % express feature points in chaser frame at the specified pose
    rMat = f_stateToPosChaserFrame(xVec, rCamVec, rFeaMat);

    %% obtain measurements
    if options.generate_measurements == 1,

        % generate measurements
        yVec = f_generateMeasurements(rMat, params.cam.f);

        % add noise to measurements
        yVec = yVec + deg2rad(0.5)*randn(size(yVec));

    else
        error('TODO: not implemented');
    end

    %% LM
    % initialise initial state estimate
    %        [x; y; z; phi; theta; psi];
    xHatVec0 = xHatMat(idx,:)';
    %xHatVec0(1:3) = [0;0;30];
    %xHatVec0(4:6) = deg2rad(0)*ones(3,1);

    if idx > 3, params.lm.reinit_att_noise_std = deg2rad(2); end
    
    % find pose estimate
    xHatVec = f_LM_adaptive_reinit(xHatVec0,yVec,rCamVec,rFeaMat,params.lm);
    
    xHatMat(idx,:) = xHatVec';
    xHatMat(idx+1,:) = xHatVec';
    
    % find conjugate pose estimate
    xHatVec2 = f_findConjPose(xHatVec);
    % refine conjugate pose estimate
    xHatVec2 = f_LM_adaptive(xHatVec2,yVec,rCamVec,rFeaMat,params.lm);
    
    rMatHat = f_stateToPosChaserFrame(xHatVec, rCamVec, rFeaMat);
    rMatHat2 = f_stateToPosChaserFrame(xHatVec2, rCamVec, rFeaMat);

    %% plot
    % plot feature point positions wrt chaser in chaser frame
    if options.plot == 1,
        fig = figure;
        
        f_posePlotter(fig, rMat, 'truth');
    %end
    %if options.plot == 1,
        f_posePlotter(fig, rMatHat, 'estimate');
        %f_posePlotter(fig, rMatHat2, 'est_conj');
    end
    
    if options.write_to_file == 1,
        fprintf(fileID_write,fmt_write,xHatVec');
        fprintf(fileID_write2,fmt_write,xHatVec2');
    end
end

if options.time == 1, toc; end

if options.write_to_file == 1,
    fclose(fileID_write);
end

post_analysis_tracking;