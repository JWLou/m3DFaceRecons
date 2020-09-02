function [J_rColor, rColor] = jacobian_rColor(params, img)
% Calculate the color residuals with respect to the current parameters and the Jacobian matrix 
% of color residuals with respect to the parameters based on IRLS solver - [alpha, beta, delta, 
% vecAngle(phi(angleX), theta(angleY), psi(angleZ)), T, f, (u0, v0), gamma(3*9)]

%% Parameters
% camera parameters
phi = params.vecAngle(1); theta = params.vecAngle(2); psi = params.vecAngle(3);
Rx = [1, 0, 0; 0, cos(phi), -sin(phi); 0, sin(phi), cos(phi)];
Ry = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
Rz = [cos(psi), -sin(psi), 0; sin(psi), cos(psi), 0; 0, 0, 1];
R = Rz*Ry*Rx; 
T = params.T;
f = params.f; u0 = params.u0; v0 = params.v0;
% shape model parameters
alpha = params.alpha; delta = params.delta; 
aid = params.mu_id+params.mu_exp; Eid = params.pc_id; Eexp = params.pc_exp;
Mgeo = aid + Eid*alpha + Eexp*delta; mPts = length(aid)/3; Mgeo = reshape(Mgeo, [3, mPts]);
matTri = params.tri(:, params.idxFaceTri); mTri = size(matTri, 2);
% reflectance(albedo) model parameters
aalb = params.mu_tex; Ealb = params.pc_tex; beta = params.beta;
Malb = aalb + Ealb*beta; Malb = reshape(Malb, [3, mPts]);
% illumination model parameters
gamma = params.gamma;
gammaR = gamma(1:9); gammaG = gamma(10:18); gammaB = gamma(19:27);
% the histogram of parameters
vec_mParam = params.vec_mParam; mParam = sum(vec_mParam);
mAlpha = vec_mParam(1); mBeta = vec_mParam(2); mDelta = vec_mParam(3); 
mR = vec_mParam(4); mT = vec_mParam(5); mIntrinCam = vec_mParam(6);
m1 = mAlpha; m2 = m1+mBeta; m3 = m2+mDelta; m4 = m3+mR; m5 = m4+mT; m6 = m5+mIntrinCam;
% image
imgW = size(img, 2); imgH = size(img, 1); 
% coefficents used in SH basis functions
c0 = 1/sqrt(4*pi); c1 = sqrt(3/(4*pi)); c2 = 3*sqrt(5/(12*pi)); 
C = [c0; c1; c1; c1; c2; c2; c2; c2/2; c2/(2*sqrt(3))];
% partial derivatives of the rotation matrix with respect to rotation angles
Rx_phi = [0, 0, 0; 0, -sin(phi), -cos(phi); 0, cos(phi), -sin(phi)];
Ry_theta = [-sin(theta), 0, cos(theta); 0, 0, 0; -cos(theta), 0, -sin(theta)];
Rz_psi = [-sin(psi), -cos(psi), 0; cos(psi), -sin(psi), 0; 0, 0, 0];
dR_dphi = Rz*Ry*Rx_phi;
dR_dtheta = Rz*Ry_theta*Rx;
dR_dpsi = Rz_psi*Ry*Rx;

%% Rasterization
% camera coordinates
MgeoC = R*Mgeo+repmat(T, [1, mPts]);
% uv coordinates
matPuv = [f*MgeoC(1, :)./MgeoC(3, :)+u0*ones(1, mPts); f*MgeoC(2, :)./MgeoC(3, :)+v0*ones(1, mPts)];
% z-buffer
zBuffer = zeros(imgH, imgW);
% synthesized image
imgSyn = zeros(size(img));
% mapping between triangles on the mesh and image pixels
mapTriPix = zeros(imgH, imgW);
baryCoord = cell(imgH, imgW); % barycentric coordinates of image pixels covered by projected 3D face triangles
for k = 1:mTri
    % traverse all the triangles on the face mesh
    idxP1 = matTri(1, k); idxP2 = matTri(2, k); idxP3 = matTri(3, k);    
    uvTri = matPuv(:, [idxP1, idxP2, idxP3]);
    uvTri = uvTri';
    
    % test if the projected triangle is reasonable and falls inside the image
    uvTriMin = ceil(min(uvTri)); uvTriMax = floor(max(uvTri)); 
    
    % calculate the triangle normal
    Pc1 = MgeoC(:, idxP1); Pc2 = MgeoC(:, idxP2); Pc3 =  MgeoC(:, idxP3);   
    U = Pc3-Pc1; V = Pc2-Pc1;    
    normTri = [U(2)*V(3)-U(3)*V(2); U(3)*V(1)-U(1)*V(3); U(1)*V(2)-U(2)*V(1)];
    unitNorm = normTri/norm(normTri); % normalise the original normal to a unit vector
    nx = unitNorm(1); ny = unitNorm(2); nz = unitNorm(3); 
    % spherical harmonics basis functions
    Y = [1; nx; ny; nz; nx*ny; nx*nz; ny*nz; nx^2-ny^2; 3*nz^2-1];
    Ynm = C.*Y; 
    
    % a visible triangle should be in front of the camera
    if normTri(3) <= 0
        % traverse image pixels fall inside the triangle
        uvTriMin = max([1, 1; uvTriMin]);
        uvTriMax = min([imgW, imgH; uvTriMax]); 
        for u = uvTriMin(1):uvTriMax(1)
            for v = uvTriMin(2):uvTriMax(2)
                [isInside, lambdaP1, lambdaP2, lambdaP3] = pointInTri([u, v], uvTri); % test if the pixel is covered by the triangle
                Xs = lambdaP1*Pc1+lambdaP2*Pc2+lambdaP1*Pc3; 
                if isInside == true && Xs(3) >= params.f
                    if mapTriPix(v, u) == 0 || zBuffer(v, u) > Xs(3)
                        zBuffer(v, u) = Xs(3);
                        mapTriPix(v, u) = k; % assign the triangle index
                        % calculate the triangle texture value
                        albTri = lambdaP1*Malb(:, idxP1)+lambdaP2*Malb(:, idxP2)+lambdaP3*Malb(:, idxP3);    
                        imgSyn(v, u, :) = [albTri(1)*gammaR'*Ynm, albTri(2)*gammaG'*Ynm, albTri(3)*gammaB'*Ynm]; % assign pixel color values
                        baryCoord{v, u} = [lambdaP1, lambdaP2, lambdaP3]; % assign barycentric coordinates
                    end      
                end
            end
        end
    end
end

%% Calculate the color residuals with respect to the current parameters and the Jacobian matrix of residual functions with respect to the color(photo) consistency term
[GuR, GvR] = imgradientxy(img(:, :, 1)); % image gradient 
[GuG, GvG] = imgradientxy(img(:, :, 2));
[GuB, GvB] = imgradientxy(img(:, :, 3));
[rows, cols] = find(mapTriPix); % visible triangles
mPixFace = length(rows); 
J_rColor = zeros(3*mPixFace, mParam);
rColor = zeros(3*mPixFace, 1);
for i = 1:mPixFace
    % find the 3D mesh triangle contributes to the image pixel
    idxRow = rows(i); idxCol = cols(i);
    idxTri = mapTriPix(idxRow, idxCol);
    idxP1 = matTri(1, idxTri); idxP2 = matTri(2, idxTri); idxP3 = matTri(3, idxTri); 
    idxP1ACs = idxP1*3-2:idxP1*3; idxP2ACs = idxP2*3-2:idxP2*3; idxP3ACs = idxP3*3-2:idxP3*3;
    Pw1 = Mgeo(:, idxP1); Pw2 = Mgeo(:, idxP2); Pw3 =  Mgeo(:, idxP3); 
    Pc1 = MgeoC(:, idxP1); Pc2 = MgeoC(:, idxP2); Pc3 =  MgeoC(:, idxP3); 
    lambdaTri = baryCoord{idxRow, idxCol}; % barycentric coordinates
    
    % calculate the triangle normal
    U = Pc3-Pc1; V = Pc2-Pc1; 
    normTri = [U(2)*V(3)-U(3)*V(2); U(3)*V(1)-U(1)*V(3); U(1)*V(2)-U(2)*V(1)]; % the original triangle normal, norm = [nx, ny, nz]
    nx = normTri(1); ny = normTri(2); nz = normTri(3); 
    lenNorm = norm(normTri); lenNorm2 = lenNorm^2;
    unitNorm = normTri/lenNorm; % normalise the original normal to a unit vector
    nxU = unitNorm(1); nyU = unitNorm(2); nzU = unitNorm(3); 
     
    % calculate the albedo
    albTri = lambdaTri*[Malb(:, idxP1)'; Malb(:, idxP2)'; Malb(:, idxP3)'];   
    cR = albTri(1); cG = albTri(2); cB = albTri(3);
    idxR = i*3-2; idxG = idxR+1; idxB = idxR+2; % indices of R, G, B channels
    
    % spherical harmonics basis functions
    Y = [1; nxU; nyU; nzU; nxU*nyU; nxU*nzU; nyU*nzU; nxU^2-nyU^2; 3*nzU^2-1];
    Ynm = C.*Y; 
      
    %% color resduials 
    rColor(idxR) = cR*gammaR'*Ynm-double(img(idxRow, idxCol, 1)); 
    rColor(idxG) = cG*gammaG'*Ynm-double(img(idxRow, idxCol, 2));
    rColor(idxB) = cB*gammaB'*Ynm-double(img(idxRow, idxCol, 3));
       
    %% the IRLS weight
    wIRLS_rColor = 1/norm([rColor(idxR), rColor(idxG), rColor(idxB)]); 
    wIRLS_rColor_sqrt = sqrt(wIRLS_rColor);
    
    %% partial derivatives with respect to alpha    
    dPw1_dalpha = Eid(idxP1ACs, :); dPc1_dalpha = R*dPw1_dalpha;
    dPw2_dalpha = Eid(idxP2ACs, :); dPc2_dalpha = R*dPw2_dalpha;
    dPw3_dalpha = Eid(idxP3ACs, :); dPc3_dalpha = R*dPw3_dalpha;

    % C_S - colour synthesized
    dU_dalpha = dPc3_dalpha-dPc1_dalpha;
    dV_dalpha = dPc2_dalpha-dPc1_dalpha;
    dnx_dalpha = dU_dalpha(2, :)*V(3)+U(2)*dV_dalpha(3, :)-dU_dalpha(3, :)*V(2)-U(3)*dV_dalpha(2, :);
    dny_dalpha = dU_dalpha(3, :)*V(1)+U(3)*dV_dalpha(1, :)-dU_dalpha(1, :)*V(3)-U(1)*dV_dalpha(3, :);
    dnz_dalpha = dU_dalpha(1, :)*V(2)+U(1)*dV_dalpha(2, :)-dU_dalpha(2, :)*V(1)-U(2)*dV_dalpha(1, :);
    dlenNorm_dalpha = (nx*dnx_dalpha+ny*dny_dalpha+nz*dnz_dalpha)/lenNorm;
    dunitNormX_dalpha = (dnx_dalpha*lenNorm-nx*dlenNorm_dalpha)/lenNorm2;
    dunitNormY_dalpha = (dny_dalpha*lenNorm-ny*dlenNorm_dalpha)/lenNorm2;
    dunitNormZ_dalpha = (dnz_dalpha*lenNorm-nz*dlenNorm_dalpha)/lenNorm2;
    dYnm_dalpha = repmat(C, [1, mAlpha]).*[zeros(1, mAlpha); dunitNormX_dalpha; dunitNormY_dalpha; dunitNormZ_dalpha; ...
                      dunitNormX_dalpha*nyU+nxU*dunitNormY_dalpha; dunitNormX_dalpha*nzU+nxU*dunitNormZ_dalpha; dunitNormY_dalpha*nzU+nyU*dunitNormZ_dalpha; ...
                      2*nxU*dunitNormX_dalpha-2*nyU*dunitNormY_dalpha; ...
                      6*nzU*dunitNormZ_dalpha];
    dCsR_dalpha = cR*gammaR'*dYnm_dalpha; % R
    dCsG_dalpha = cG*gammaG'*dYnm_dalpha; % G
    dCsB_dalpha = cB*gammaB'*dYnm_dalpha; % B
     
    % C_I - colour input
    dPuv1_dalpha = (f/Pc1(3)^2)*[dPc1_dalpha(1, :)*Pc1(3)-Pc1(1)*dPc1_dalpha(3, :); ...
                                 dPc1_dalpha(2, :)*Pc1(3)-Pc1(2)*dPc1_dalpha(3, :)];
    dPuv2_dalpha = (f/Pc2(3)^2)*[dPc2_dalpha(1, :)*Pc2(3)-Pc2(1)*dPc2_dalpha(3, :); ...
                                 dPc2_dalpha(2, :)*Pc2(3)-Pc2(2)*dPc2_dalpha(3, :)];
    dPuv3_dalpha = (f/Pc3(3)^2)*[dPc3_dalpha(1, :)*Pc3(3)-Pc3(1)*dPc3_dalpha(3, :); ...
                                 dPc3_dalpha(2, :)*Pc3(3)-Pc3(2)*dPc3_dalpha(3, :)];
    dcenTriu_dalpha = lambdaTri*[dPuv1_dalpha(1, :); dPuv2_dalpha(1, :); dPuv3_dalpha(1, :)];  
    dcenTriv_dalpha = lambdaTri*[dPuv1_dalpha(2, :); dPuv2_dalpha(2, :); dPuv3_dalpha(2, :)]; 
    dCiR_dalpha = GuR(idxRow, idxCol)*dcenTriu_dalpha+GvR(idxRow, idxCol)*dcenTriv_dalpha; % R
    dCiG_dalpha = GuG(idxRow, idxCol)*dcenTriu_dalpha+GvG(idxRow, idxCol)*dcenTriv_dalpha; % G
    dCiB_dalpha = GuB(idxRow, idxCol)*dcenTriu_dalpha+GvB(idxRow, idxCol)*dcenTriv_dalpha; % B

    % fill in the Jacobian matrix
    J_rColor(idxR, 1:m1) = dCsR_dalpha-dCiR_dalpha; % R
    J_rColor(idxG, 1:m1) = dCsG_dalpha-dCiG_dalpha; % G
    J_rColor(idxB, 1:m1) = dCsB_dalpha-dCiB_dalpha; % B
    
    %% partial derivatives with respect to beta
    Ealb_P1 = Ealb(idxP1ACs, :);
    Ealb_P2 = Ealb(idxP2ACs, :);
    Ealb_P3 = Ealb(idxP3ACs, :);
    dCsR_dbeta = lambdaTri*[Ealb_P1(1, :); Ealb_P2(1, :); Ealb_P3(1, :)];
    dCsG_dbeta = lambdaTri*[Ealb_P1(2, :); Ealb_P2(2, :); Ealb_P3(2, :)];
    dCsB_dbeta = lambdaTri*[Ealb_P1(3, :); Ealb_P2(3, :); Ealb_P3(3, :)];
    
    % fill in the Jacobian matrix
    J_rColor(idxR, m1+1:m2) = gammaR'*Ynm*dCsR_dbeta; % R
    J_rColor(idxG, m1+1:m2) = gammaG'*Ynm*dCsG_dbeta; % G
    J_rColor(idxB, m1+1:m2) = gammaB'*Ynm*dCsB_dbeta; % B
    
    %% partial derivatives with respect to delta 
    dPw1_ddelta = Eexp(idxP1ACs, :); dPc1_ddelta = R*dPw1_ddelta;
    dPw2_ddelta = Eexp(idxP2ACs, :); dPc2_ddelta = R*dPw2_ddelta;
    dPw3_ddelta = Eexp(idxP3ACs, :); dPc3_ddelta = R*dPw3_ddelta;

    % C_S - colour synthesized
    dU_ddelta = dPc3_ddelta-dPc1_ddelta;
    dV_ddelta = dPc2_ddelta-dPc1_ddelta;
    dnx_ddelta = dU_ddelta(2, :)*V(3)+U(2)*dV_ddelta(3, :)-dU_ddelta(3, :)*V(2)-U(3)*dV_ddelta(2, :);
    dny_ddelta = dU_ddelta(3, :)*V(1)+U(3)*dV_ddelta(1, :)-dU_ddelta(1, :)*V(3)-U(1)*dV_ddelta(3, :);
    dnz_ddelta = dU_ddelta(1, :)*V(2)+U(1)*dV_ddelta(2, :)-dU_ddelta(2, :)*V(1)-U(2)*dV_ddelta(1, :);
    dlenNorm_ddelta = (nx*dnx_ddelta+ny*dny_ddelta+nz*dnz_ddelta)/lenNorm;
    dunitNormX_ddelta = (dnx_ddelta*lenNorm-nx*dlenNorm_ddelta)/lenNorm2;
    dunitNormY_ddelta = (dny_ddelta*lenNorm-ny*dlenNorm_ddelta)/lenNorm2;
    dunitNormZ_ddelta = (dnz_ddelta*lenNorm-nz*dlenNorm_ddelta)/lenNorm2;
    dYnm_ddelta = repmat(C, [1, mDelta]).*[zeros(1, mDelta); dunitNormX_ddelta; dunitNormY_ddelta; dunitNormZ_ddelta; ...
                      dunitNormX_ddelta*nyU+nxU*dunitNormY_ddelta; dunitNormX_ddelta*nzU+nxU*dunitNormZ_ddelta; dunitNormY_ddelta*nzU+nyU*dunitNormZ_ddelta; ...
                      2*nxU*dunitNormX_ddelta-2*nyU*dunitNormY_ddelta; ...
                      6*nzU*dunitNormZ_ddelta];
    dCsR_ddelta = cR*gammaR'*dYnm_ddelta; % R
    dCsG_ddelta = cG*gammaG'*dYnm_ddelta; % G
    dCsB_ddelta = cB*gammaB'*dYnm_ddelta; % B
     
    % C_I - colour input
    dPuv1_ddelta = (f/Pc1(3)^2)*[dPc1_ddelta(1, :)*Pc1(3)-Pc1(1)*dPc1_ddelta(3, :); ...
                                 dPc1_ddelta(2, :)*Pc1(3)-Pc1(2)*dPc1_ddelta(3, :)];
    dPuv2_ddelta = (f/Pc2(3)^2)*[dPc2_ddelta(1, :)*Pc2(3)-Pc2(1)*dPc2_ddelta(3, :); ...
                                 dPc2_ddelta(2, :)*Pc2(3)-Pc2(2)*dPc2_ddelta(3, :)];
    dPuv3_ddelta = (f/Pc3(3)^2)*[dPc3_ddelta(1, :)*Pc3(3)-Pc3(1)*dPc3_ddelta(3, :); ...
                                 dPc3_ddelta(2, :)*Pc3(3)-Pc3(2)*dPc3_ddelta(3, :)];
    dcenTriu_ddelta = lambdaTri*[dPuv1_ddelta(1, :); dPuv2_ddelta(1, :); dPuv3_ddelta(1, :)];  
    dcenTriv_ddelta = lambdaTri*[dPuv1_ddelta(2, :); dPuv2_ddelta(2, :); dPuv3_ddelta(2, :)]; 
    dCiR_ddelta = GuR(idxRow, idxCol)*dcenTriu_ddelta+GvR(idxRow, idxCol)*dcenTriv_ddelta; % R
    dCiG_ddelta = GuG(idxRow, idxCol)*dcenTriu_ddelta+GvG(idxRow, idxCol)*dcenTriv_ddelta; % G
    dCiB_ddelta = GuB(idxRow, idxCol)*dcenTriu_ddelta+GvB(idxRow, idxCol)*dcenTriv_ddelta; % B

    % fill in the Jacobian matrix
    J_rColor(idxR, m2+1:m3) = dCsR_ddelta-dCiR_ddelta; % R
    J_rColor(idxG, m2+1:m3) = dCsG_ddelta-dCiG_ddelta; % G
    J_rColor(idxB, m2+1:m3) = dCsB_ddelta-dCiB_ddelta; % B

    %% partial derivatives with respect to rotation angles
    dPc1_dangle = [dR_dphi*Pw1, dR_dtheta*Pw1, dR_dpsi*Pw1];
    dPc2_dangle = [dR_dphi*Pw2, dR_dtheta*Pw2, dR_dpsi*Pw2];
    dPc3_dangle = [dR_dphi*Pw3, dR_dtheta*Pw3, dR_dpsi*Pw3];
    
    % C_S - colour synthesized
    dU_dangle = dPc3_dangle-dPc1_dangle;
    dV_dangle = dPc2_dangle-dPc1_dangle;   
    dnx_dangle = dU_dangle(2, :)*V(3)+U(2)*dV_dangle(3, :)-dU_dangle(3, :)*V(2)-U(3)*dV_dangle(2, :);
    dny_dangle = dU_dangle(3, :)*V(1)+U(3)*dV_dangle(1, :)-dU_dangle(1, :)*V(3)-U(1)*dV_dangle(3, :);
    dnz_dangle = dU_dangle(1, :)*V(2)+U(1)*dV_dangle(2, :)-dU_dangle(2, :)*V(1)-U(2)*dV_dangle(1, :);
    dlenNorm_dangle = (nx*dnx_dangle+ny*dny_dangle+nz*dnz_dangle)/lenNorm;
    dunitNormX_dangle = (dnx_dangle*lenNorm-nx*dlenNorm_dangle)/lenNorm2;
    dunitNormY_dangle = (dny_dangle*lenNorm-ny*dlenNorm_dangle)/lenNorm2;
    dunitNormZ_dangle = (dnz_dangle*lenNorm-nz*dlenNorm_dangle)/lenNorm2;
    dYnm_dangle = repmat(C, [1, mR]).*[zeros(1, mR); dunitNormX_dangle; dunitNormY_dangle; dunitNormZ_dangle; ...
                      dunitNormX_dangle*nyU+nxU*dunitNormY_dangle; dunitNormX_dangle*nzU+nxU*dunitNormZ_dangle; dunitNormY_dangle*nzU+nyU*dunitNormZ_dangle; ...
                      2*nxU*dunitNormX_dangle-2*nyU*dunitNormY_dangle; ...
                      6*nzU*dunitNormZ_dangle];
    dCsR_dangle = cR*gammaR'*dYnm_dangle; % R
    dCsG_dangle = cG*gammaG'*dYnm_dangle; % G
    dCsB_dangle = cB*gammaB'*dYnm_dangle; % B
    
    % C_I - colour input    
    dPuv1_dangle = (f/Pc1(3)^2)*[dPc1_dangle(1, :)*Pc1(3)-Pc1(1)*dPc1_dangle(3, :); ...
                                 dPc1_dangle(2, :)*Pc1(3)-Pc1(2)*dPc1_dangle(3, :)];
    dPuv2_dangle = (f/Pc2(3)^2)*[dPc2_dangle(1, :)*Pc2(3)-Pc2(1)*dPc2_dangle(3, :); ...
                                 dPc2_dangle(2, :)*Pc2(3)-Pc2(2)*dPc2_dangle(3, :)];
    dPuv3_dangle = (f/Pc3(3)^2)*[dPc3_dangle(1, :)*Pc3(3)-Pc3(1)*dPc3_dangle(3, :); ...
                                 dPc3_dangle(2, :)*Pc3(3)-Pc3(2)*dPc3_dangle(3, :)];
    dcenTriu_dangle = lambdaTri*[dPuv1_dangle(1, :); dPuv2_dangle(1, :); dPuv3_dangle(1, :)];  
    dcenTriv_dangle = lambdaTri*[dPuv1_dangle(2, :); dPuv2_dangle(2, :); dPuv3_dangle(2, :)]; 
    dCiR_dangle = GuR(idxRow, idxCol)*dcenTriu_dangle+GvR(idxRow, idxCol)*dcenTriv_dangle; % R
    dCiG_dangle = GuG(idxRow, idxCol)*dcenTriu_dangle+GvG(idxRow, idxCol)*dcenTriv_dangle; % G
    dCiB_dangle = GuB(idxRow, idxCol)*dcenTriu_dangle+GvB(idxRow, idxCol)*dcenTriv_dangle; % B
    
    % fill in the Jacobian matrix
    J_rColor(idxR, m3+1:m4) = dCsR_dangle-dCiR_dangle; % R
    J_rColor(idxG, m3+1:m4) = dCsG_dangle-dCiG_dangle; % G
    J_rColor(idxB, m3+1:m4) = dCsB_dangle-dCiB_dangle; % B

    %% partial derivatives with respect to the translation
    dPuv1_dT = [f/Pc1(3), 0, -f*Pc1(1)/(Pc1(3)^2); 0, f/Pc1(3), -f*Pc1(2)/(Pc1(3)^2)]; 
    dPuv2_dT = [f/Pc2(3), 0, -f*Pc2(1)/(Pc2(3)^2); 0, f/Pc2(3), -f*Pc2(2)/(Pc2(3)^2)];
    dPuv3_dT = [f/Pc3(3), 0, -f*Pc3(1)/(Pc3(3)^2); 0, f/Pc3(3), -f*Pc3(2)/(Pc3(3)^2)];   
    dcenTriu_dT = lambdaTri*[dPuv1_dT(1, :); dPuv2_dT(1, :); dPuv3_dT(1, :)];
    dcenTriv_dT = lambdaTri*[dPuv1_dT(2, :); dPuv2_dT(2, :); dPuv3_dT(2, :)];
    dCiR_dT = GuR(idxRow, idxCol)*dcenTriu_dT+GvR(idxRow, idxCol)*dcenTriv_dT; % R
    dCiG_dT = GuG(idxRow, idxCol)*dcenTriu_dT+GvG(idxRow, idxCol)*dcenTriv_dT; % G
    dCiB_dT = GuB(idxRow, idxCol)*dcenTriu_dT+GvB(idxRow, idxCol)*dcenTriv_dT; % B
    
    % fill in the Jacobian matrix
    J_rColor(idxR, m4+1:m5) = -dCiR_dT; % R
    J_rColor(idxG, m4+1:m5) = -dCiG_dT; % G
    J_rColor(idxB, m4+1:m5) = -dCiB_dT; % B

    %% partial derivatives with respect to projection parameters
    dPuv1_dproj = [Pc1(1)/Pc1(3), 1, 0; Pc1(2)/Pc1(3), 0, 1];
    dPuv2_dproj = [Pc2(1)/Pc2(3), 1, 0; Pc2(2)/Pc2(3), 0, 1];
    dPuv3_dproj = [Pc3(1)/Pc3(3), 1, 0; Pc3(2)/Pc3(3), 0, 1];
    dcenTriu_dproj = lambdaTri*[dPuv1_dproj(1, :); dPuv2_dproj(1, :); dPuv3_dproj(1, :)];
    dcenTriv_dproj = lambdaTri*[dPuv1_dproj(2, :); dPuv2_dproj(2, :); dPuv3_dproj(2, :)];
    dCiR_dproj = GuR(idxRow, idxCol)*dcenTriu_dproj+GvR(idxRow, idxCol)*dcenTriv_dproj; % R
    dCiG_dproj = GuG(idxRow, idxCol)*dcenTriu_dproj+GvG(idxRow, idxCol)*dcenTriv_dproj; % G
    dCiB_dproj = GuB(idxRow, idxCol)*dcenTriu_dproj+GvB(idxRow, idxCol)*dcenTriv_dproj; % B 
    
    % fill the Jacobian matrix
    J_rColor(idxR, m5+1:m6) = -dCiR_dproj; % R
    J_rColor(idxG, m5+1:m6) = -dCiG_dproj; % G
    J_rColor(idxB, m5+1:m6) = -dCiB_dproj; % B
    
    %% partial derivatives with respect to gamma
    J_rColor(idxR, m6+1:m6+9) = cR*Ynm';   % R
    J_rColor(idxG, m6+10:m6+18) = cG*Ynm'; % G
    J_rColor(idxB, m6+19:m6+27) = cB*Ynm'; % B
    
    %% assign the IRLS weight
    J_rColor(idxR:idxB, :) = wIRLS_rColor_sqrt*J_rColor(idxR:idxB, :); 
    rColor(idxR:idxB) = wIRLS_rColor_sqrt*rColor(idxR:idxB);
end

% assign the normalization factor
wColor_sqrt_norm = sqrt(1/mPixFace);
J_rColor = wColor_sqrt_norm*J_rColor;
rColor = wColor_sqrt_norm*rColor;

end