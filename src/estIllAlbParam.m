function params = estIllAlbParam(params, img)
% Estimate initial illumination and albedo parameters

%% Parameters
% camera parameters
phi = params.vecAngle(1); theta = params.vecAngle(2); psi = params.vecAngle(3);
R = rotationMatrix(phi, theta, psi); T = params.T; f = params.f; u0 = params.u0; v0 = params.v0;
% shape model parameters
alpha = params.alpha; delta = params.delta; 
aid = params.mu_id+params.mu_exp; Eid = params.pc_id; Eexp = params.pc_exp;
Mgeo = aid + Eid*alpha + Eexp*delta; mPts = length(aid)/3; Mgeo = reshape(Mgeo, [3, mPts]);
matTri = params.tri; mTri = size(matTri, 2);
% reflectance(albedo) model parameters
aalb = params.mu_tex; Ealb = params.pc_tex; beta = params.beta;
Malb = aalb+Ealb*beta; Malb = reshape(Malb, [3, mPts]); 
mBeta = params.vec_mParam(2); 
% image
imgW = size(img, 2); imgH = size(img, 1); 
% coefficents used in SH basis functions
c0 = 1/sqrt(4*pi); c1 = sqrt(3/(4*pi)); c2 = 3*sqrt(5/(12*pi)); 
C = [c0; c1; c1; c1; c2; c2; c2; c2/2; c2/(2*sqrt(3))];

%% Rasterization
% camera coordinates
MgeoC = R*Mgeo+repmat(T, [1, mPts]);
% uv coordinates
matPuv = [f*MgeoC(1, :)./MgeoC(3, :)+u0*ones(1, mPts); f*MgeoC(2, :)./MgeoC(3, :)+v0*ones(1, mPts)];
% z-buffer
zBuffer = zeros(imgH, imgW);
% mapping between triangles on the mesh and image pixels
mapTriPix = zeros(imgH, imgW);
baryCoord = cell(imgH, imgW); % barycentric coordinates of image pixels covered by projected 3D face triangles
for k = 1:mTri
    % traverse all the triangles on the face mesh
    idxP1 = matTri(1, k); idxP2 = matTri(2, k); idxP3 = matTri(3, k);    
    uvTri = matPuv(:, [idxP1, idxP2, idxP3]);
    uvTri = uvTri';
    
    % test if the projected triangle intersects with the image
    uvTriMin = ceil(min(uvTri)); uvTriMax = floor(max(uvTri)); 
    if uvTriMax(1) < uvTriMin(1) || uvTriMax(2) < uvTriMin(2) || uvTriMax(1) > imgW || uvTriMin(1) < 0 || uvTriMax(2) > imgH || uvTriMin(2) < 0
        continue;
    end
    
    % the center of the triangle 
    Pc1 = MgeoC(:, idxP1); Pc2 = MgeoC(:, idxP2); Pc3 =  MgeoC(:, idxP3);   
    U = Pc3-Pc1; V = Pc2-Pc1;    
    % calculate the triangle normal
    normTri = [U(2)*V(3)-U(3)*V(2); U(3)*V(1)-U(1)*V(3); U(1)*V(2)-U(2)*V(1)];
    
    % a visible triangle should be in front of the camera
    if normTri(3) <= 0
        % traverse image pixels fall inside the triangle
        uvTriMin = max([1, 1; uvTriMin]);
        uvTriMax = min([imgW, imgH; uvTriMax]); 
        for u = uvTriMin(1):uvTriMax(1)
            for v = uvTriMin(2):uvTriMax(2)
                [isInside, lambdaP1, lambdaP2, lambdaP3] = pointInTri([u, v], uvTri); % test if the pixel is covered by the triangle
                if isInside == true
                    Xs = lambdaP1*Pc1+lambdaP2*Pc2+lambdaP1*Pc3; 
                    if (mapTriPix(v, u) == 0 || zBuffer(v, u) > Xs(3))
                        zBuffer(v, u) = Xs(3);
                        mapTriPix(v, u) = k; % assign the triangle index
                        baryCoord{v, u} = [lambdaP1, lambdaP2, lambdaP3]; % assign barycentric coordinates
                    end      
                end
            end
        end
    end
end

%% Calculate the color residuals with respect to the current parameters and the Jacobian matrix of residual functions with respect to the color(photo) consistency term
[rows, cols] = find(mapTriPix); % visible triangles
mPixFace = length(rows); 
AgammaR = zeros(mPixFace, 9); AgammaG = zeros(mPixFace, 9); AgammaB = zeros(mPixFace, 9); 
bgammaR = zeros(mPixFace, 1); bgammaG = zeros(mPixFace, 1); bgammaB = zeros(mPixFace, 1);
AalbR = zeros(mPixFace, mBeta); AalbG = zeros(mPixFace, mBeta); AalbB = zeros(mPixFace, mBeta); 
balbR = zeros(mPixFace, 1); balbG = zeros(mPixFace, 1); balbB = zeros(mPixFace, 1);
matY = zeros(9, mPixFace); 
for i = 1:mPixFace
    % find the 3D mesh triangle contributes to the image pixel
    idxRow = rows(i);
    idxCol = cols(i);
    idxTri = mapTriPix(idxRow, idxCol);
    idxP1 = matTri(1, idxTri); idxP2 = matTri(2, idxTri); idxP3 = matTri(3, idxTri);  
    Pc1 = MgeoC(:, idxP1); Pc2 = MgeoC(:, idxP2); Pc3 =  MgeoC(:, idxP3); 
         
    % calculate the albedo
    lambdaTri = baryCoord{idxRow, idxCol}; % barycentric coordinates
    albTri = lambdaTri*[Malb(:, idxP1)'; Malb(:, idxP2)'; Malb(:, idxP3)'];   
    cR = albTri(1); cG = albTri(2); cB = albTri(3);
    
    % spherical harmonics basis functions
    U = Pc3-Pc1; V = Pc2-Pc1; 
    normTri = [U(2)*V(3)-U(3)*V(2); U(3)*V(1)-U(1)*V(3); U(1)*V(2)-U(2)*V(1)]; % the original triangle normal, norm = [nx, ny, nz]
    unitNorm = normTri/norm(normTri); % normalise the original normal to a unit vector
    nxU = unitNorm(1); nyU = unitNorm(2); nzU = unitNorm(3); 
    Y = [1; nxU; nyU; nzU; nxU*nyU; nxU*nzU; nyU*nzU; nxU^2-nyU^2; 3*nzU^2-1]; matY(:, i) = Y;
    Ynm = C.*Y; 
    
    % fitting targets
    bgammaR(i) = double(img(idxRow, idxCol, 1)); 
    bgammaG(i) = double(img(idxRow, idxCol, 2));
    bgammaB(i) = double(img(idxRow, idxCol, 3));
    balbR(i) = bgammaR(i); 
    balbG(i) = bgammaG(i); 
    balbB(i) = bgammaB(i);
        
    % data points
    AgammaR(i, :) = cR*Ynm'; % R
    AgammaG(i, :) = cG*Ynm'; % G
    AgammaB(i, :) = cB*Ynm'; % B
    
    idxRP1 = idxP1*3-2; idxGP1 = idxRP1+1; idxBP1 = idxRP1+2; 
    idxRP2 = idxP2*3-2; idxGP2 = idxRP2+1; idxBP2 = idxRP2+2; 
    idxRP3 = idxP3*3-2; idxGP3 = idxRP3+1; idxBP3 = idxRP3+2; 
    AalbR(i, :) = lambdaTri*Ealb([idxRP1, idxRP2, idxRP3], :);
    AalbG(i, :) = lambdaTri*Ealb([idxGP1, idxGP2, idxGP3], :); 
    AalbB(i, :) = lambdaTri*Ealb([idxBP1, idxBP2, idxBP3], :);
end

% estimate illumination parameters
gammaR = (AgammaR'*AgammaR)\AgammaR'*bgammaR;
gammaG = (AgammaG'*AgammaG)\AgammaG'*bgammaG;
gammaB = (AgammaB'*AgammaB)\AgammaB'*bgammaB;
params.gamma(1:9) = gammaR';
params.gamma(10:18) = gammaG';
params.gamma(19:27) = gammaB';

% estimate albedo parameters
Aalb = [repmat((gammaR'*matY)', [1, mBeta]).*AalbR; repmat((gammaG'*matY)', [1, mBeta]).*AalbG; repmat((gammaB'*matY)', [1, mBeta]).*AalbB];
balb = [balbR; balbG; balbB];
params.beta = (Aalb'*Aalb)\Aalb'*balb;

end
