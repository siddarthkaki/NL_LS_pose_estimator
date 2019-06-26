function [xHatVec] = f_lsqnonlin_reinit(funhandle,xHatVec0,params)
%F_LSQNONLIN_REINIT Multiple Random Reinitialisation Wrapper for MATLAB's
%                   in-built lsqnonlin solver

xHatMat = zeros(params.lm.num_init,6);
JVec = zeros(params.lm.num_init,1);

for init_idx = 1:params.lm.num_init,
    
    xHatVec0new = xHatVec0;
    xHatVec0new(4:6) = normrnd(0, params.lm.reinit_att_noise_std, [3 1]);
    
    [xHatVec, Jval] = lsqnonlin(funhandle,xHatVec0new,[],[],params.lsqnonlin);
   
    % store estimate and corresponding cost for ith reinitialisation
    xHatMat(init_idx,:) = xHatVec';
    JVec(init_idx) = Jval;
end

% return estimate with the lowest corresponding cost
[~,min_J_idx] = min(JVec);
xHatVec = xHatMat(min_J_idx,:)';

end