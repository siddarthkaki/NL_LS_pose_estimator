function [yVec] = f_generateMeasurements(rMat,focal_length)
%F_GENERATEMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here

[numPts,~] = size(rMat);

% project feature points to image plane
imgPtMat = zeros(numPts,2);
for idx = 1:numPts,
    rVeci = rMat(idx,:)';
    imgPti = f_3Dto2D_projection(rVeci,focal_length);
    imgPtMat(idx,:) = imgPti';
end

% azimuth & elevation measurements for feature points
azVec = zeros(numPts,1);
elVec = zeros(numPts,1);
for idx = 1:numPts,
    imgPti = imgPtMat(idx,:)';
    
    azVec(idx) = atan2(imgPti(1),focal_length);
    elVec(idx) = atan2(imgPti(2),focal_length);
end

yVec = []; % vector of measurements
for idx = 1:numPts,
    yVec = [yVec; azVec(idx); elVec(idx)];
end

end

