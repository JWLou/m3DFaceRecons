function visualizeNormals
% Visualize vertex normals of on the 3D face 

%% Parameters
% shape model parameters
modelExp = load('./model/modelExp.mat');
modelId = load('./model/modelId.mat');
Mgeo = double(modelId.mu_id + modelExp.mu_exp);
mPts = length(Mgeo)/3;
Mgeo = reshape(Mgeo, [3, mPts]);
matTri = modelId.tri; 
mTri = size(matTri, 2);

%% Draw
trisurf(matTri', Mgeo(1, :), Mgeo(2, :), Mgeo(3, :), 0, 'edgecolor', 'none');
light('Position', [0 0 10^8], 'Style', 'infinite');
lighting gouraud
trimesh(matTri', Mgeo(1, :), Mgeo(2, :), Mgeo(3, :));
axis equal
axis vis3d
rotate3d on;
title('3D shape');
xlabel('x');
ylabel('y');
zlabel('z');
grid off;
axis on;
view([0 90]);
hold on;
for i = 1:mTri
    idxP1 = matTri(1, i); idxP2 = matTri(2, i); idxP3 = matTri(3, i);    
    Pw1 = Mgeo(:, idxP1); Pw2 = Mgeo(:, idxP2); Pw3 =  Mgeo(:, idxP3); 
    cTri = double((Pw1+Pw2+Pw3)/3);
    
    % calculate and draw the normal of each triangle
    U = Pw3-Pw1; V = Pw2-Pw1; 
    normTri = [U(2)*V(3)-U(3)*V(2); U(3)*V(1)-U(1)*V(3); U(1)*V(2)-U(2)*V(1)];
    normTri = 5*10^4*normTri/norm(normTri); % scale the norm to unified length
    quiver3(cTri(1), cTri(2), cTri(3), normTri(1), normTri(2), normTri(3), 0.5);
end
hold off;

end