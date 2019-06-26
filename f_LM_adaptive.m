function [xHatVec] = f_LM_adaptive(xHatVec0,yVec,rCamVec,rFeaMat,params)
%F_LM_ADAPTIVE Non-linear Least-Squares Levenberg–Marquardt Solver
%              for Pose with Adaptive Marquardt Parameter Selection
%
% Measurements: Relative bearing (azimuth, elevation) to feature points with 
%               known correspondences 
%
% Outputs:      Relative pose between chaser frame and target frame
%               (x, y, z, roll, pitch, yaw)

% LM iteration
xHatVec = xHatVec0;

lambda = params.lambda;

count = 0;
while true,
    
    xHatVec_old = xHatVec;

    H    = f_HMat(xHatVec, rCamVec, rFeaMat);
    
    [~, Hy] = size(H);
    
    M = H.'*H;

    % check for positive definiteness of H.'*H
    try chol(M);
        %disp('Matrix is symmetric positive definite.')
    catch
        lambda = 1;%norm(H)*0.1;
        %disp('Matrix is not symmetric positive definite')
    end
    
    yHatVec = f_yHatVec(xHatVec, rCamVec, rFeaMat);
    Jold = f_J(yVec,yHatVec);
    
    %dxHatVec = (M + lambda*eye(Hy))\H.'*(yVec - yHatVec);
    
    while true,
        
        %dxHatVec = (H.'*H + lambda*eye(Hy))\H.'*(yVec - yHatVec);
        dxHatVec = (M + lambda*diag(diag(M)))\H.'*(yVec - yHatVec);

        Jnew = f_J(yVec,f_yHatVec(xHatVec + dxHatVec, rCamVec, rFeaMat));

        if Jnew >= Jold,
            lambda = lambda*5;%max(2*lambda,norm(H)*0.1);
        else
            lambda = lambda/5;
            break;
        end
    end
    
    xHatVec = xHatVec + dxHatVec;
    
    count = count + 1;
    
    if count >= params.max_count,
        break;
    end
    
    if norm(xHatVec - xHatVec_old) <= params.eps,
        break;
    end
end
end

%% Dependent functions
% Jacobian matrix
function HMat = f_HMat(xHatVec, rCamVec, rFeaMat)

[numPts,~] = size(rFeaMat); 

rRelVec = xHatVec(1:3);
phi = xHatVec(4);
theta = xHatVec(5);
psi = xHatVec(6);

HMat = [];

for idx = 1:numPts,
    rFeaVeci = rFeaMat(idx,:)';
    HVecsi = [ (rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))/((rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2 + (rRelVec(1) - rCamVec(1) + rFeaVeci(1)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - rFeaVeci(3)*cos(phi)*sin(theta))^2),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             0, -(rRelVec(1) - rCamVec(1) + rFeaVeci(1)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - rFeaVeci(3)*cos(phi)*sin(theta))/((rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2 + (rRelVec(1) - rCamVec(1) + rFeaVeci(1)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - rFeaVeci(3)*cos(phi)*sin(theta))^2),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ((rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi) - rFeaVeci(1)*cos(phi)*sin(psi))*(rFeaVeci(1)*cos(psi) - rCamVec(1)*cos(theta) + rRelVec(1)*cos(theta) + rFeaVeci(2)*sin(psi) - rCamVec(3)*sin(theta) + rRelVec(3)*sin(theta)))/((rRelVec(3) - rCamVec(3) + rFeaVeci(3)*cos(phi)*cos(theta) + rFeaVeci(1)*cos(psi)*sin(theta) + rFeaVeci(2)*sin(psi)*sin(theta) - rFeaVeci(2)*cos(psi)*cos(theta)*sin(phi) + rFeaVeci(1)*cos(theta)*sin(phi)*sin(psi))^2 + (rRelVec(1) - rCamVec(1) + rFeaVeci(1)*cos(psi)*cos(theta) - rFeaVeci(3)*cos(phi)*sin(theta) + rFeaVeci(2)*cos(theta)*sin(psi) + rFeaVeci(2)*cos(psi)*sin(phi)*sin(theta) - rFeaVeci(1)*sin(phi)*sin(psi)*sin(theta))^2), -(((rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))/(rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta)) + ((rFeaVeci(1)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - rFeaVeci(3)*cos(phi)*sin(theta))*(rRelVec(1) - rCamVec(1) + rFeaVeci(1)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - rFeaVeci(3)*cos(phi)*sin(theta)))/(- rCamVec(3) + rRelVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2)*(rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2)/((rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2 + (rRelVec(1) - rCamVec(1) + rFeaVeci(1)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - rFeaVeci(3)*cos(phi)*sin(theta))^2), -(((rFeaVeci(1)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - rFeaVeci(2)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)))/(rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta)) - ((rFeaVeci(1)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) - rFeaVeci(2)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)))*(rRelVec(1) - rCamVec(1) + rFeaVeci(1)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - rFeaVeci(3)*cos(phi)*sin(theta)))/(- rCamVec(3) + rRelVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2)*(rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2)/((rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2 + (rRelVec(1) - rCamVec(1) + rFeaVeci(1)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - rFeaVeci(3)*cos(phi)*sin(theta))^2);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0, (rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))/((rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2 + (rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi) - rFeaVeci(1)*cos(phi)*sin(psi))^2),                                                                                                                                                                   -(rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi) - rFeaVeci(1)*cos(phi)*sin(psi))/((rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2 + (rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi) - rFeaVeci(1)*cos(phi)*sin(psi))^2), (((rFeaVeci(3)*cos(phi) - rFeaVeci(2)*cos(psi)*sin(phi) + rFeaVeci(1)*sin(phi)*sin(psi))/(rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta)) + ((rFeaVeci(3)*cos(theta)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi)*cos(theta) - rFeaVeci(1)*cos(phi)*cos(theta)*sin(psi))*(rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi) - rFeaVeci(1)*cos(phi)*sin(psi)))/(- rCamVec(3) + rRelVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2)*(rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2)/((rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2 + (rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi) - rFeaVeci(1)*cos(phi)*sin(psi))^2),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           -((rFeaVeci(1)*(cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) - rFeaVeci(3)*cos(phi)*sin(theta))*(rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi) - rFeaVeci(1)*cos(phi)*sin(psi)))/((rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2 + (rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi) - rFeaVeci(1)*cos(phi)*sin(psi))^2),                                                                                                                                                                                                                                         -(((rFeaVeci(1)*cos(phi)*cos(psi) + rFeaVeci(2)*cos(phi)*sin(psi))/(rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta)) - ((rFeaVeci(1)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) - rFeaVeci(2)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)))*(rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi) - rFeaVeci(1)*cos(phi)*sin(psi)))/(- rCamVec(3) + rRelVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2)*(rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2)/((rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(phi)*cos(theta))^2 + (rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(phi)*cos(psi) - rFeaVeci(1)*cos(phi)*sin(psi))^2) ];
    HMat = [HMat; HVecsi];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
end
end

% measurement estimates
function yHatVec = f_yHatVec(xHatVec, rCamVec, rFeaMat)

[numPts,~] = size(rFeaMat); 

rRelVec = xHatVec(1:3);
phi = xHatVec(4);
theta = xHatVec(5);
psi = xHatVec(6);

yHatVec = zeros(numPts*2,1);
for idx = 1:numPts,
    rFeaVeci = rFeaMat(idx,:)';
    yHatVec(2*idx-1) = atan2(rRelVec(1) - rCamVec(1) + rFeaVeci(1)*(cos(psi)*cos(theta) - sin(psi)*sin(theta)*sin(phi)) + rFeaVeci(2)*(cos(theta)*sin(psi) + cos(psi)*sin(theta)*sin(phi)) - rFeaVeci(3)*sin(theta)*cos(phi), rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(psi)*sin(phi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(theta)*cos(phi));
    yHatVec(2*idx)   =                                                                                  atan2(rRelVec(2) - rCamVec(2) + rFeaVeci(3)*sin(phi) + rFeaVeci(2)*cos(psi)*cos(phi) - rFeaVeci(1)*sin(psi)*cos(phi), rRelVec(3) - rCamVec(3) + rFeaVeci(1)*(cos(psi)*sin(theta) + cos(theta)*sin(psi)*sin(phi)) + rFeaVeci(2)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(3)*cos(theta)*cos(phi));
end
end

% compute LS cost
function J = f_J(yVec,yHatVec)
    J = (yVec-yHatVec).'*(yVec-yHatVec);
end