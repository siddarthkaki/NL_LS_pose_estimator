delete data/pose_trajectory_true.txt
fileID = fopen('data/pose_trajectory_true.txt','a');
fmt = '%d %d %d %d %d %d\n';

pos = [0;0;20];
rot = zeros(3,1);

for idx = 1:100,
    xVec = [pos; rot];
    fprintf(fileID,fmt,xVec');
    
    pos = pos + [0.05; -0.05; 0.1];
    rot = rot + deg2rad(10)*ones(3,1);
end
fclose(fileID);