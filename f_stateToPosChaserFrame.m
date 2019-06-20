function rMat_c = f_stateToPosChaserFrame(stateVec, rCamVec, rFeaMat)
%F_STATETOPOSCHASERFRAME convert state (truth or estimate) to feature point 
%                        positions wrt chaser in chaser frame

rRelVec = stateVec(1:3);
phi = stateVec(4);
theta = stateVec(5);
psi = stateVec(6);

% DCM from target frame to chaser frame
DCM = euler2dcm([phi theta psi]); % 3-1-2 sequence

% number of feature points
[numPts,~] = size(rFeaMat);

% position vectors of feature points wrt chaser in chaser frame
rMat_c = zeros(numPts,3);
for idx = 1:numPts,
    rFeaVeci = rFeaMat(idx,:)';
    rFeaVeci_c = DCM*rFeaVeci;
    rVeci = rRelVec - rCamVec + rFeaVeci_c; % position vector of feature point i wrt chaser in chaser frame
    rMat_c(idx,:) = rVeci';
end
end

