%% Relative Pose Estimator Script
% Author:   Siddarth Kaki

%% housekeeping
clear; close all; clc;

%% options
options.plot = 0;
options.generate_measurements = 1;
options.write_to_file = 1;
options.time = 1;
%options.lsqnonlin = optimoptions('lsqnonlin','Display','iter');
options.lsqnonlin.Algorithm = 'levenberg-marquardt';
options.lsqnonlin.FunctionTolerance = 1e-15;
options.jacobian_type = 'approx';
%options.jacobian_type = 'true';
%#ok<*NOCOL>
   
%% params
% measurement params
params.meas.std = deg2rad(1);

% camera params
params.cam.f = 5.5*10^(-3); % camera focal length

% LM params
params.lm.lambda = 5;
params.lm.max_count = 10000;
params.lm.max_adaptive_count = 5;
params.lm.eps = 0.01;
params.lm.num_init = 10;
params.lm.reinit_att_noise_std = 2;
params.lm.meas.std = params.meas.std;

%% init
rng(1)

%% construct feature points
% specify rigid position vector of feature points wrt target in target frame
rFeaMat = [ 0, 0, 0.5;
            0, 0,-1.5;
            0, 1, 1;
            0,-1, 1; ];
        
% rFeaMat = [ 0, 0, 0;
%             0, 0,-1;
%             0, 1.5, 1;
%             0,-1.5, 1; ];

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

filename_read = 'data/poses_true.txt';
xMat = read_poses(filename_read,1,linecount(filename_read)); % read in multiple poses to run

%xHatMat = zeros(size(xMat));

if options.write_to_file == 1,
    % setup output file to write to
    delete data/poses_estimated.txt
    fileID_write = fopen('data/poses_estimated.txt','a');
    
    delete data/conj_poses_estimated.txt
    fileID_write2 = fopen('data/conj_poses_estimated.txt','a');
    fmt_write = '%d %d %d %d %d %d\n';
end

if options.time == 1, tic; end

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
        yVec = yVec + normrnd(0,params.meas.std,size(yVec));

    else
        error('TODO: not implemented');
    end

    %% LM
    % initialise initial state estimate
    %        [x; y; z; phi; theta; psi];
    xNoise = [normrnd(0,2,[2 1]); normrnd(0,3,[1 1]); normrnd(0,pi,[3 1])];
    xHatVec0 = xVec + xNoise;
    
    % true Jacobian (slower but more accurate)
    if strcmp(options.jacobian_type, 'true'),
        
        % find pose estimate
        %xHatVec = f_LM_adaptive_reinit(xHatVec0,yVec,rCamVec,rFeaMat,params.lm);
        xHatVec = f_LM_sqrt_adaptive_reinit(xHatVec0,yVec,rCamVec,rFeaMat,params.lm);
        
        % find conjugate pose estimate
        xHatVec2 = f_findConjPose(xHatVec);
        
        % refine conjugate pose estimate
        xHatVec2 = f_LM_adaptive(xHatVec2,yVec,rCamVec,rFeaMat,params.lm);
        
    % approx Jacobian (faster but less accurate)
    elseif strcmp(options.jacobian_type, 'approx'),
        
        f_measResidWrapper = @(xHatVecParam) f_measResid(xHatVecParam, yVec, rCamVec, rFeaMat);
        
        % find pose estimate
        xHatVec = lsqnonlin(f_measResidWrapper,xHatVec0,[],[],options.lsqnonlin);
        
        % find conjugate pose estimate
        xHatVec2 = f_findConjPose(xHatVec);
        
        % refine conjugate pose estimate
        xHatVec2 = lsqnonlin(f_measResidWrapper,xHatVec2,[],[],options.lsqnonlin);
    end
    
    rMatHat = f_stateToPosChaserFrame(xHatVec, rCamVec, rFeaMat);
    rMatHat2 = f_stateToPosChaserFrame(xHatVec2, rCamVec, rFeaMat);

    error = xVec - xHatVec;
    error_deg = [error(1:3); rad2deg(error(4:6))]';

    %% plot
    % plot feature point positions wrt chaser in chaser frame
    if options.plot == 1,
        fig = figure;
        
        f_posePlotter(fig, rMat, 'truth');
    %end
    %if options.plot == 1,
        f_posePlotter(fig, rMatHat, 'estimate');
        f_posePlotter(fig, rMatHat2, 'est_conj');
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

post_analysis;