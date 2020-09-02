function drawHead(vertex, tri, keyPoints)

trisurf(tri', vertex(1, :), vertex(2, :), vertex(3, :), 0, 'edgecolor', 'none');

if nargin == 3
    keyPoints = double(keyPoints);
    hold on
    plot3(keyPoints(1, :), keyPoints(2, :), keyPoints(3, :), 'b.', 'MarkerSize', 17);    
%     msg = cell(1, size(keyPoints, 2));
%     for i = 1 : size(keyPoints, 2)
%         msg{i} = sprintf('%d', i);
%     end
%     text(keyPoints(1, :), keyPoints(2, :), keyPoints(3, :) + 5e4, msg);
end

light('Position', [0 0 -10^0], 'Style', 'infinite');
lighting gouraud
axis equal
axis vis3d
rotate3d on;
xlabel('x');
ylabel('y');
zlabel('z');
grid off;
axis on;
view([0 -90]);

% hold on;
% mTri = size(tri, 2);
% for i = 1:mTri
%     idxP1 = tri(1, i); idxP2 = tri(2, i); idxP3 = tri(3, i);    
%     Pw1 = vertex(:, idxP1); Pw2 = vertex(:, idxP2); Pw3 =  vertex(:, idxP3); 
%     cTri = (Pw1+Pw2+Pw3)/3;
%     
%     % calculate the normal
%     U = Pw3-Pw1; V = Pw2-Pw1; 
%     normTri = [U(2)*V(3)-U(3)*V(2); U(3)*V(1)-U(1)*V(3); U(1)*V(2)-U(2)*V(1)];
%     normTri = 1*10^4*normTri/norm(normTri); % scale the norm to unified length
%     % draw the normal of each triangle
%     quiver3(cTri(1), cTri(2), cTri(3), normTri(1), normTri(2), normTri(3), 0.5);
% end
% hold off;

% Rotate the 3D shape in a loop
% for i=1:72
%     camorbit(5,0,'camera');
%     M=getframe(gcf);
%     nn=frame2im(M);
%     [nn,cm]=rgb2ind(nn,256);
%     if i==1
%         imwrite(nn,cm,'out.gif','gif','LoopCount',inf,'DelayTime',0.1);
%     else
%         imwrite(nn,cm,'out.gif','gif','WriteMode','append','DelayTime',0.1);
%     end
% end

end