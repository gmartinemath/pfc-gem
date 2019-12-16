%Returns the statistics associated to each grain
%Gabriel Martine
%Last updated: Apr 12 2019

function stats = GetGrainStats(file)
	load(file);
	invalid_grains = [];
	
	%Simple statistics
	area = g_area;
	perim = g_perim;
	ipr = IPR(g_area, g_perim);
	
	%Convex hull ratio
	if exist('g_hull_ratio', 'var')
		chr = g_hull_ratio;
	else
		%Use the actual average d for the atomic distribution by assigning equal to hexagonal area to all atoms
		chr = ((g_hull_area + sqrt(2*L^2/(sqrt(3)*n_true))/2 * g_hull_perim) - g_area)./g_area;
	end
	
	%Coordination number
	coord = sum(g_adjmat, 2);
	
	%Average of neighbors' coordination number
	avg_neighbor_coord = 0*coord;
	for cc = 1:length(coord)
		avg_neighbor_coord(cc) = mean(coord(g_adjmat(cc, :)));
	end
		
	%Remove invalid and boundary grains
	invalid_grains = [find(area == 0); find(perim == 0); find(chr < 0)];
	if exist('g_on_boundary', 'var')
		invalid_grains = unique([invalid_grains; find(g_on_boundary)]);
	else
		invalid_grains = unique(invalid_grains);
	end
	area(invalid_grains) = []; perim(invalid_grains) = [];
	ipr(invalid_grains) = []; chr(invalid_grains) = []; coord(invalid_grains) = [];
	avg_neighbor_coord(invalid_grains) = [];
	
	%Package into a structure
	stats.area = area; stats.perim = perim; stats.ipr = ipr; stats.chr = chr; stats.coord = coord;
	stats.avg_neighbor_coord = avg_neighbor_coord;
end
