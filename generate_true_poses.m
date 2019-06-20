delete data/poses_true.txt
fileID = fopen('data/poses_true.txt','a');
fmt = '%d %d %d %d %d %d\n';

% pos = [0;0;25] + [randn(2,1)*1; randn(1,1)*3];
% rot = deg2rad(0)*zeros(3,1) + randn(3,1)*deg2rad(30);

for idx = 1:500,
    pos = [0;0;25] + [randn(2,1)*1; randn(1,1)*3];
    rot = deg2rad(0)*zeros(3,1) + randn(3,1)*deg2rad(90);
    xVec = [pos; rot];
    fprintf(fileID,fmt,xVec');
end
fclose(fileID);