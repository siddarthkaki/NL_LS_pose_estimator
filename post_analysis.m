%% Post Analysis

%% housekeeping
%clear; close all; clc;

%% init
filename_read_true = 'data/poses_true.txt';
xMat = read_poses(filename_read_true,1,linecount(filename_read_true));

filename_read_estimated = 'data/poses_estimated.txt';
xHatMat = read_poses(filename_read_estimated,1,linecount(filename_read_estimated));

filename_read_conj_estimated = 'data/conj_poses_estimated.txt';
xHatMat2 = read_poses(filename_read_conj_estimated,1,linecount(filename_read_conj_estimated));

numPoses = size(xMat,1);

%%
ErrMat_rad = xMat - xHatMat;
ErrMat_deg = ErrMat_rad;
ErrMat_deg(:,4:6) = rad2deg(ErrMat_rad(:,4:6));

pos_score = 100*ones(numPoses,1);
att_score = 100*ones(numPoses,1);

for idx = 1:numPoses,
    
    % pos score
    rVec = xMat(idx,1:3);
    rHatVec = xHatMat(idx,1:3);
    
    pos_score(idx) = norm(rVec - rHatVec)/norm(rVec);
    
    % att score
    tru_quat = dcm2quat(euler2dcm(xMat(idx,4:6)));
    
    est_dcm  = euler2dcm(xHatMat(idx,4:6));
    est_dcm2 = euler2dcm(xHatMat2(idx,4:6));
    
    est_quat  = dcm2quat(euler2dcm(xHatMat(idx,4:6)));
    est_quat2 = dcm2quat(euler2dcm(xHatMat2(idx,4:6)));

    tru_quat = tru_quat/quatnorm(tru_quat);
    est_quat = est_quat/quatnorm(est_quat);
    est_quat2 = est_quat2/quatnorm(est_quat2);
    
    z_quat = quatmultiply(est_quat,quatconj(tru_quat));
    z_quat2 = quatmultiply(est_quat2,quatconj(tru_quat));
    
    f_att_score = @(quat_w) 2*acosd(norm(quat_w));
    
    att_score(idx) = f_att_score(z_quat(1));
    att_score(idx) = min(f_att_score(z_quat(1)), f_att_score(z_quat2(1)));
end

% STDErr_rad = std(ErrMat_rad);
% 
MeanErr_deg = mean(ErrMat_deg);
% STDErr_deg = std(ErrMat_deg);

pos_score_mean = mean(pos_score);
pos_score_std = std(pos_score);

att_score_mean = mean(att_score);
att_score_std = std(att_score);

fprintf('\npos_score - mean: %.5f\t|\tstd: %.3f [m]\n',pos_score_mean,pos_score_std)
fprintf('\natt_score - mean: %.5f\t|\tstd: %.3f [deg]\n',att_score_mean,att_score_std)
