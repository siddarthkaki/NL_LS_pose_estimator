function [] = f_posePlotter(fig_handle, rMat_c, type)
%F_POSEPLOTTER plot feature points wrt chaser in chaser frame for visualisation
%   Detailed explanation goes here

figure(fig_handle);

title('truth=b, estimate=r')
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');

[numPts,~] = size(rMat_c);

if strcmp(type,'truth'),
    for idx = 1:numPts,
        plot3(subs(rMat_c(idx,1)), subs(rMat_c(idx,2)), subs(rMat_c(idx,3)),'*b')
        hold on
    end

    for idx = 2:numPts,
        plot3([rMat_c(1,1), rMat_c(idx,1)],[rMat_c(1,2), rMat_c(idx,2)],[rMat_c(1,3), rMat_c(idx,3)],'b')
        hold on
    end

elseif strcmp(type,'estimate'),
    for idx = 1:numPts,
        plot3(subs(rMat_c(idx,1)), subs(rMat_c(idx,2)), subs(rMat_c(idx,3)),'*r')
        hold on
    end

    for idx = 2:numPts,
        plot3([rMat_c(1,1), rMat_c(idx,1)],[rMat_c(1,2), rMat_c(idx,2)],[rMat_c(1,3), rMat_c(idx,3)],'r')
        hold on
    end
    
elseif strcmp(type,'est_conj'),
    title('truth=b, estimate=r, conj\_est=g')
    for idx = 1:numPts,
        plot3(subs(rMat_c(idx,1)), subs(rMat_c(idx,2)), subs(rMat_c(idx,3)),'*g')
        hold on
    end

    for idx = 2:numPts,
        plot3([rMat_c(1,1), rMat_c(idx,1)],[rMat_c(1,2), rMat_c(idx,2)],[rMat_c(1,3), rMat_c(idx,3)],'g')
        hold on
    end 
end
grid on
axis equal
end