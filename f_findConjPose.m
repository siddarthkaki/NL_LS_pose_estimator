function [xHatVec2] = f_findConjPose(xHatVec)
%F_FINDCONJPOSE Summary of this function goes here
%   Detailed explanation goes here
    est_dcm = euler2dcm(xHatVec(4:6));

    p0 = est_dcm*[0;0;0];
    p1 = est_dcm*[0;0;1];
    p2 = est_dcm*[0;1;0];

    p01 = p1-p0;
    p02 = p2-p0;
    n = cross(p01,p02);
    if n(3) < 0, n = -n; end

    phi  = acos(dot(n,[0;0;1])/norm(n));
    aHat = cross(n,[0;0;1]); aHat = aHat/norm(aHat);
    est_dcm2 = rotationMatrix(aHat,-2*phi)*est_dcm;
    
    xHatVec2 = xHatVec;
    xHatVec2(4:6) = dcm2euler(est_dcm2);

end