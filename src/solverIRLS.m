function params = solverIRLS(params, img, lan2d)
% Apply Iteratively Reweighted Least Squares(IRLS) solver to minimize the objective 
% function and update the parameters - [alpha, beta, delta, vecAngle(phi(angleX),
% theta(angleY), psi(angleZ)), T, f, (u0, v0), gamma(3*9)]

fprintf('The IRLS step ...\n');

% term weights
wColor = params.wColor; wColor_sqrt = sqrt(wColor);
wLan = params.wLan; wLan_sqrt = sqrt(wLan);
wReg = params.wReg; wReg_sqrt = sqrt(wReg);

% the histogram of parameters
vec_mParam = params.vec_mParam; mParam = sum(vec_mParam);
mAlpha = vec_mParam(1); mBeta = vec_mParam(2); mDelta = vec_mParam(3); 
mR = vec_mParam(4); mT = vec_mParam(5); mIntrinCam = vec_mParam(6);
m1 = mAlpha; m2 = m1+mBeta; m3 = m2+mDelta; m4 = m3+mR; m5 = m4+mT; m6 = m5+mIntrinCam;

% show initial results before the Gauss-Newton process
% showRes(params, img, lan2d);
for i = 1:params.nIterIRLS     
    %% Gauss-Newton step
     [J_rColor, rColor] = jacobian_rColor(params, img); % photo consistency 
     [J_rLan, rLan] = jacobian_rLan(params, lan2d); % feature alignment 
     [J_rReg, rReg] = jacobian_rReg(params); % statistical regularization
     
     % print errors
     fprintf('The %dth Gauss-Newton step: \n', i-1);
     fprintf('Color residuals: %.4f \n', norm(rColor));
     fprintf('Landmark residuals: %.4f \n', norm(rLan));
     fprintf('Regularization residuals: %.4f \n', norm(rReg));
     fprintf('----------------------------- \n');
     
     % assign the energy term weight
     J_rColor = wColor_sqrt*J_rColor; rColor = wColor_sqrt*rColor;
     J_rLan = wLan_sqrt*J_rLan; rLan = wLan_sqrt*rLan;
     J_rReg = wReg_sqrt*J_rReg; rReg = wReg_sqrt*rReg;
     
     % the overall Jacobian and residuals
     J_rAll = [J_rColor; J_rLan; J_rReg]; 
     rAll = [rColor; rLan; rReg];
     
%      J_rAll(:, m5+1:m5+3) = [];
%      deltaParam = zeros(mParam-3, 1); 
%      params.gamma = params.gamma+deltaParam(m5+1:end);
     
     % call the PCG solover to compute the update of parameters 
     A = J_rAll'*J_rAll; b = -J_rAll'*rAll; 
     deltaParam = zeros(mParam, 1); 
     deltaParam = solverPCG(A, b, deltaParam, params.nIterPCG);
     
     % update parameters
     params.alpha = params.alpha+deltaParam(1:m1);
     params.beta = params.beta+deltaParam(m1+1:m2);
     params.delta = params.delta+deltaParam(m2+1:m3);
     params.vecAngle = params.vecAngle+deltaParam(m3+1:m4);
     params.T = params.T+deltaParam(m4+1:m5);
     params.f = params.f+deltaParam(m5+1);
     params.u0 = params.u0+deltaParam(m5+2);
     params.v0 = params.v0+deltaParam(m5+3);
     params.gamma = params.gamma+deltaParam(m6+1:end);

     % show updated results after one iteration 
%      showRes(params, img, lan2d); 
end

end