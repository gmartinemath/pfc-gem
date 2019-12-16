%This computes area and interface lengths by summing the area of the relevant atom's Voronoi regions
%Grain area is simply the sum of the areas of the atoms in one grain
%Grain perimeter is half of the sum of the areas of the atoms at the grain boundary and those that touch the grain
%For the interface lengths, atoms at triple and higher junctions count in all three (or more) interfaces
%Gabriel Martine
%Last updated: April 25 2019

function [g_area, g_perim, g_adjmat, int_length, mis_angle] = MeasureGrains(g_angle, ...
	atoms_grain_id, ghost_ind, n_true, grain_n, adj_matrix, atoms_area, interatomic_distance)
	
	%Initialize blank return arrays
	g_adjmat = false(grain_n);
	g_area = zeros(grain_n, 1);
	g_perim = g_area;
	int_length = zeros(grain_n);
	
	for cur_grain = 1:grain_n
		%Fetch all true atoms that are part of the current grain and compute its area
		fetch_id = find(atoms_grain_id(1:n_true) == cur_grain);
		g_area(cur_grain) = sum(atoms_area(fetch_id));
		
		%Initialize the array of atoms that have already been used to compute perimeters
		atoms_done = false(n_true, 1);
		
		%Now check if the grain's atoms touch a different grain
		for check_i = fetch_id
			atom_neighbors = ghost_ind(adj_matrix{check_i});
			
			%Remove all neighbors that are part of the current grain
			atom_neighbors = atom_neighbors(atoms_grain_id(atom_neighbors) ~= cur_grain);
			
			if ~isempty(atom_neighbors)
				%If the current atom touches a neighboring grain, add its Voronoi area to the perimeter
				g_perim(cur_grain) = g_perim(cur_grain) + atoms_area(check_i);
				
				%Go over the neighbor grains to add neighbor areas at the correct place
				for neighbor_grain = unique(atoms_grain_id(atom_neighbors))
					g_adjmat(cur_grain, neighbor_grain) = true;
					int_length(cur_grain, neighbor_grain) = int_length(cur_grain, neighbor_grain) + atoms_area(check_i);
					
					%Add the contribution of the neighbors if they haven't been considered yet
					for neighbor_atom = atom_neighbors
						if ~atoms_done(neighbor_atom)
							g_perim(cur_grain) = g_perim(cur_grain) + atoms_area(neighbor_atom);
							int_length(cur_grain, neighbor_grain)=int_length(cur_grain, neighbor_grain)+atoms_area(neighbor_atom);
							atoms_done(neighbor_atom) = true;
						end
					end
				end
			end
		end
	end
	
	%Recover the perimeter from the measure of the boundary area
	g_perim = g_perim / (2*interatomic_distance);
	int_length = int_length / (2*interatomic_distance);
	
	%Compute the angle and area discrepancy matrices
	mis_angle = zeros(grain_n);
	[rr, cc] = find(g_adjmat);
	for ind = 1:length(rr)
		mis_angle(rr(ind),cc(ind)) = g_angle(rr(ind)) - g_angle(cc(ind));
	end
	mis_angle = AngleDist(mis_angle, pi/3);
end
