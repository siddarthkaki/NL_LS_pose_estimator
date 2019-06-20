function [point2DVec] = f_3Dto2D_projection(point3DVec,f)
%F_3DTO2D_PROJECTION Summary of this function goes here
%   Detailed explanation goes here

PMat = [f 0 0 0; 0 f 0 0; 0 0 1 0];

projVec = PMat*[point3DVec; 1];

x = projVec(1)/projVec(3);
y = projVec(2)/projVec(3);

if isnan(x), x = 0; end
if isnan(y), y = 0; end

point2DVec = [x;y];

end

