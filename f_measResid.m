function [measResidVec] = f_measResid(xHatVec, yVec, rCamVec, rFeaMat)
%F_ Summary of this function goes here
%   Detailed explanation goes here

[numPts,~] = size(rFeaMat); 

rRelVec = xHatVec(1:3);
phi = xHatVec(4);
theta = xHatVec(5);
psi = xHatVec(6);

yHatVec = zeros(numPts*2,1);
for yHat_idx = 1:numPts,
    rFeaVeci = rFeaMat(yHat_idx,:)';
    yHatVec(2*yHat_idx-1) = atan2(rRelVec(1) - rCamVec(1) + rFeaVeci(1)*(cos(psi)*cos(theta) - sin(psi)*sin(theta)*sin(phi)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(theta)*sin(phi)) - rFeaVeci(3)*sin(theta)*cos(phi), rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(psi)*sin(phi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(theta)*cos(phi));
    yHatVec(2*yHat_idx)   =                                                                                  atan2(rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(psi)*cos(phi) - rFeaVeci(1)*sin(psi)*cos(phi), rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(psi)*sin(phi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(theta)*cos(phi));
end

measResidVec = yVec - yHatVec;

end


