function idxFaceTri = extractTri()

%% Parameters
% shape model parameters
modelExp = load('./model/modelExp.mat');
modelId = load('./model/modelId.mat');
symlist_tri = load('./model/symlist_tri.mat');
Mgeo = double(modelId.mu_id + modelExp.mu_exp);
mPts = length(Mgeo)/3;
Mgeo = reshape(Mgeo, [3, mPts]);
matTri = modelId.tri; 
keyPoints = Mgeo(:, modelId.idxKpts);

%% Draw
trimesh(matTri(:, symlist_tri(:, nonfaceTri))', Mgeo(1, :), Mgeo(2, :), Mgeo(3, :));
keyPoints = double(keyPoints);
hold on
plot3(keyPoints(1, :), keyPoints(2, :), keyPoints(3, :), 'b.', 'MarkerSize', 17); 
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
% length(symlist_tri(2, :))
for rang = 52886:52937
rang
idxTri = symlist_tri(2, rang);
idxPts = matTri(:, idxTri);
msg = sprintf('%d', rang);
cenTri = mean(Mgeo(:, idxPts'), 2);
text(cenTri(1), cenTri(2), cenTri(3) + 5e2, msg);
% plot3(cenTri(1), cenTri(2), cenTri(3), 'r.', 'MarkerSize', 17); 
% trimesh(idxPts', Mgeo(1, :), Mgeo(2, :), Mgeo(3, :));
pause;
end
hold off;

nonfaceTri = [16089:25361, 25465:25600, 25728:25845, 25980:26096, 26240:26332, 26507:26590, 26761:26840, ...
              42362:42483, 42578:42700, 42795:42916, 43012:43132, 43229:43347, 43443:43560, 43653:43771, ...
              43865:43980, 44075:44188, 44285:44396, 44492:44603, 44696:44809, 44901:45011, 45106:45212, ...
              45306:45411, 45503:45608, 45698:45803, 45892:45996, 46085:46187, 46276:46376, 46466:46563, ...
              46652:46748, 46837:46931, 47021:47112, 47201:47291, 47379:47468, 47556:47643, 47731:47816, ...
              47906:47987, 48076:48156, 48246:48323, 48411:48488, 48576:48651, 48736:48812, 48896:48971, ...
              49056:49128, 49211:49283, 49366:49436, 49521:49587, 49671:49736, 49821:49883, 49966:50028, ...
              50111:50171, 50256:50312, 50396:50451, 50536:50588, 50672:50723, 50806:50856, 50938:50987, ...
              51066:51117, 51196:51247, 51326:51377, 51456:51507, 51586:51637, 51716:51767, 51846:51897, ...
              51976:52027, 52106:52157, 52236:52287, 52366:52417, 52496:52547, 52626:52677, 52756:52807, ...
              52886:52937];
faceTri = setdiff([1:52937], nonfaceTri);
          
idxFaceTri = symlist_tri(:, faceTri);

end
