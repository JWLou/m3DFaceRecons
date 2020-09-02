function [imgSyn, mapTriPix] = synFaceImg(params, img, lan2d)
% Synthesize face image using the estimated parametric face model

%% Parameters
% camera parameters
phi = params.vecAngle(1); theta = params.vecAngle(2); psi = params.vecAngle(3);
R = rotationMatrix(phi, theta, psi); T = params.T;
f = params.f; u0 = params.u0; v0 = params.v0;
% shape model parameters
alpha = params.alpha; delta = params.delta; 
aid = params.mu_id+params.mu_exp; Eid = params.pc_id; Eexp = params.pc_exp;
Mgeo = aid + Eid*alpha + Eexp*delta; mPts = length(aid)/3; Mgeo = reshape(Mgeo, [3, mPts]);
matTri = params.tri; mTri = size(matTri, 2);
% reflectance(albedo) model parameters
aalb = params.mu_tex; Ealb = params.pc_tex; beta = params.beta;
Malb = aalb+Ealb*beta; Malb = reshape(Malb, [3, mPts]);
% illumination model parameters
gamma = params.gamma;
% image
imgW = size(img, 2); 
imgH = size(img, 1); 
% face patch
faceC = mean(lan2d);
faceRange = range(lan2d);
faceW = faceRange(1); faceH = faceRange(2);
faceLen = min(faceW, faceH);

%% Rasterization
% coefficents used in SH basis functions
c0 = 1/sqrt(4*pi); c1 = sqrt(3/(4*pi)); c2 = 3*sqrt(5/(12*pi)); 
C = [c0; c1; c1; c1; c2; c2; c2; c2/2; c2/(2*sqrt(3))];
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
for k = 1:mTri
    % traverse all the triangles on the face mesh
    idxP1 = matTri(1, k); idxP2 = matTri(2, k); idxP3 = matTri(3, k);    
    Pc1 = MgeoC(:, idxP1); Pc2 = MgeoC(:, idxP2); Pc3 =  MgeoC(:, idxP3);   
    cTriC = (Pc1+Pc2+Pc3)/3;
    uvTri = matPuv(:, [idxP1, idxP2, idxP3]);
    uvTri = uvTri';
    
    % test if the projected triangle falls inside the image
    uvTriMin = ceil(min(uvTri)); uvTriMax = floor(max(uvTri)); 
    if uvTriMax(1) < uvTriMin(1) || uvTriMax(2) < uvTriMin(2) || uvTriMax(1) > imgW || uvTriMin(1) < 0 || uvTriMax(2) > imgH || uvTriMin(2) < 0
        continue;
    end
    
    % calculate the triangle normal
    U = Pc3-Pc1; V = Pc2-Pc1; 
    normTri = [U(2)*V(3)-U(3)*V(2); U(3)*V(1)-U(1)*V(3); U(1)*V(2)-U(2)*V(1)];
    unitNorm = normTri/norm(normTri); % normalise the original normal to a unit vector
    nx = unitNorm(1); ny = unitNorm(2); nz = unitNorm(3); 
     
    % calculate the albedo
    albTri = (Malb(:, idxP1)+Malb(:, idxP2)+Malb(:, idxP3))/3;   
    cR = albTri(1); cG = albTri(2); cB = albTri(3);
    
    % spherical harmonics basis functions
    Y = [1; nx; ny; nz; nx*ny; nx*nz; ny*nz; nx^2-ny^2; 3*nz^2-1];
    Ynm = C.*Y; 
    gammaR = gamma(1:9); gammaG = gamma(10:18); gammaB = gamma(19:27);
    
    % calculate the triangle texture value
    tR = cR*gammaR'*Ynm; % R
    tG = cG*gammaG'*Ynm; % G
    tB = cB*gammaB'*Ynm; % B 
    
    % traverse image pixels fall inside the triangle
    uvTriMin = max([1, 1; uvTriMin]);
    uvTriMax = min([imgW, imgH; uvTriMax]); 
    for u = uvTriMin(1):uvTriMax(1)
        for v = uvTriMin(2):uvTriMax(2)
            isInside = pointInTri([u, v], uvTri); % test if the pixel is covered by the triangle
            if isInside == true
                if (imgSyn(v, u) == 0 || zBuffer(v, u) > cTriC(3))
                    zBuffer(v, u) = cTriC(3);
                    imgSyn(v, u, :) = [tR, tG, tB]; % assign pixel color values        
                    mapTriPix(v, u) = k; % assign the triangle index
                end      
            end
        end
    end
end

end