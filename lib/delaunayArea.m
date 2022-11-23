function Area = delaunayArea(DT)
% Finds the total area of a Delaunay triangulation.
%
% Input
% DT = Delaunary triangulation
%
% Output
% Area = Area of the triangulation
%%
N_triangles = size(DT.ConnectivityList,1);
Areas = nan(N_triangles,1);
Triangles_points = DT.Points(DT.ConnectivityList,:);
for i=1:N_triangles
    % Form the two vectors that defines the triangle and calculate their cross product.
    Triangle_points = Triangles_points(3*(i-1) + (1:3),:);
    Triangle_vecs = bsxfun(@minus,Triangle_points,Triangle_points(3,:));
    Areas(i) = abs(Triangle_vecs(1,1)*Triangle_vecs(2,2) - Triangle_vecs(1,2)*Triangle_vecs(2,1));
end
Area = sum(Areas);