function matVertTex = calVertTex(params)
% Calculate the texture of each vertex on the mesh

%% Parameters
% camera parameters
phi = params.vecAngle(1); theta = params.vecAngle(2); psi = params.vecAngle(3);
R = rotationMatrix(phi, theta, psi); T = params.T;
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

%% Rasterization
% coefficents used in SH basis functions
c0 = 1/sqrt(4*pi); c1 = sqrt(3/(4*pi)); c2 = 3*sqrt(5/(12*pi)); 
C = [c0; c1; c1; c1; c2; c2; c2; c2/2; c2/(2*sqrt(3))];
% camera coordinates
MgeoC = R*Mgeo+repmat(T, [1, mPts]);
% vertex normals
matNormVert = zeros(3, mPts);
countNormVert = zeros(mPts, 1);
% vertex texture
matVertTex = zeros(3, mPts);
for k = 1:mTri
    % traverse all the triangles on the face mesh
    idxP1 = matTri(1, k); idxP2 = matTri(2, k); idxP3 = matTri(3, k);    
    Pc1 = MgeoC(:, idxP1); Pc2 = MgeoC(:, idxP2); Pc3 =  MgeoC(:, idxP3);   
        
    % calculate the triangle normal
    U = Pc3-Pc1; V = Pc2-Pc1; 
    normTri = [U(2)*V(3)-U(3)*V(2); U(3)*V(1)-U(1)*V(3); U(1)*V(2)-U(2)*V(1)];
    unitNorm = normTri/norm(normTri); % normalise the original normal to a unit vector
     
    % accumulate the vertex normal
    countNormVert([idxP1, idxP2, idxP3]) =  countNormVert([idxP1, idxP2, idxP3]) + 1;
    matNormVert(:, [idxP1, idxP2, idxP3]) = matNormVert(:, [idxP1, idxP2, idxP3]) + repmat(unitNorm, [1, 3]);
end

for k = 1:mPts
    matNormVert(:, k) = matNormVert(:, k)/countNormVert(k);
    matNormVert(:, k) = matNormVert(:, k)/norm(matNormVert(:, k));
    nx = matNormVert(1, k); ny = matNormVert(2, k); nz = matNormVert(3, k); 
    % calculate the albedo 
    cR = Malb(1, k); cG = Malb(2, k); cB = Malb(3, k);
    
    % spherical harmonics basis functions
    Y = [1; nx; ny; nz; nx*ny; nx*nz; ny*nz; nx^2-ny^2; 3*nz^2-1];
    Ynm = C.*Y; 
    gammaR = gamma(1:9); gammaG = gamma(10:18); gammaB = gamma(19:27);
    
    % calculate the vertex texture value
    tR = cR*gammaR'*Ynm; % R
    tG = cG*gammaG'*Ynm; % G
    tB = cB*gammaB'*Ynm; % B 
    matVertTex(:, k) = [tR; tG; tB]; 
end

end



