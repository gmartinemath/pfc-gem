%This function takes in the (x,y) vectors of points and returns their adjacency "matrix". The matrix is
%actually a cell whose {ith} entry is the list of neighbors' ID. This is much faster than a sparse array
%Gabriel Martine
%Last updated: May 4 2017

function adj = PointAdjMatrix(x, y)
	N = length(x);
	
	%Compute the Delaunay triangulation of the points
	tri = delaunayn([x', y']);
	
	%Grow each cell entry by matching the triangulation's vertices together
	adj = cell(N, 1);
	for i = 1:size(tri, 1)
		adj{tri(i,1)} = [adj{tri(i,1)}, tri(i,2), tri(i,3)];
		adj{tri(i,2)} = [adj{tri(i,2)}, tri(i,1), tri(i,3)];
		adj{tri(i,3)} = [adj{tri(i,3)}, tri(i,1), tri(i,2)];
	end
	
	%The repetition is stripped out
	for i = 1:N
		sorted = sort(adj{i});
		adj{i} = [sorted(1), sorted(1+find(sorted(2:end) > sorted(1:end-1)))];
	end
end
