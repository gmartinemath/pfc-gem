%This function measures the convex hull area and perimeter of the region defined by atomic positions
%This does NOT loop across boundaries. Grains that cross boundaries are simply ignored for simplicity.
%Gabriel Martine
%Last updated: June 6 2018

function [g_hull_area, g_hull_perim] = MeasureConvexHull(atoms_grain_id, atoms_x, atoms_y, n_true, grain_n)
	g_hull_area = zeros(grain_n, 1);
	g_hull_perim = zeros(grain_n, 1);
	
	for cur_grain = 1:grain_n
		%Fetch all atoms that are part of the current grain and discard any grain with ghost atoms
		fetch_id = find(atoms_grain_id == cur_grain);
		if any(fetch_id > n_true)
			continue;
		end
		
		%Otherwise, compute the convex hull of the atomic positions
		g_x = atoms_x(fetch_id); g_y = atoms_y(fetch_id);
		K = convhull(g_x, g_y); g_x = g_x(K)'; g_y = g_y(K)';
		
		%Compute the hull's area and perimeter so we can approximate the grain's convex hull by area+perim*d/2
		sg_x = [g_x(end); g_x(1:end-1)]; sg_y = [g_y(end); g_y(1:end-1)];
		g_hull_area(cur_grain) = 0.5*abs(sum(g_x .* sg_y - g_y .* sg_x));
		g_hull_perim(cur_grain) = sum(sqrt((g_x-sg_x).^2 + (g_y-sg_y).^2));
	end
end

