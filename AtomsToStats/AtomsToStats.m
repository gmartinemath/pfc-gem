%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Atoms to Stats
%%%%%%%%%%%%%% This code takes in atomic positions and outputs grain statistics assuming
%%%%%%%%%%%%%% the underlying lattice is hexagonal. The method is based on finding the
%%%%%%%%%%%%%% spatial location of all atoms, resulting domain's Voronoi and Delaunay
%%%%%%%%%%%%%% tesselations and grain properties based on their atoms' local environment.
%%%%%%%%%%%%%% The domain is assumed to be a square for simplicity.
%%%%%%%%%%%%%% 
%%%%%%%%%%%%%% CLI INPUT VARIABLES
%%%%%%%%%%%%%% in_dir:		The input folder (this acts on AtomicPos_*.mat files)
%%%%%%%%%%%%%% ff_id:		The ID of the files to be converted
%%%%%%%%%%%%%% L:			The square's side
%%%%%%%%%%%%%% out_dir:		The output folder
%%%%%%%%%%%%%% ipr_factor:	The blow up factor for the Voronoi cell's isoperimetric ratio (1/gamma in the paper)
%%%%%%%%%%%%%% angle_th:	The smallest misorientation (in rad) between neighbor grains
%%%%%%%%%%%%%% grain_th:	The minimum number of atoms per identified grain (pre-filling)
%%%%%%%%%%%%%% RMN1:		True if the program should remove grains with only one neighbor
%%%%%%%%%%%%%% RMANG:		True if the program should remove grains too close in orientation to a neighbor
%%%%%%%%%%%%%% rmang_th:	The threshold angle to be used to remove neighbors too close in orientation
%%%%%%%%%%%%%% PFIG:		True if the program should save the angle and grain pictures
%%%%%%%%%%%%%% PFIG_N:		Number of pixels per dimension in the output angle and grain pictures
%%%%%%%%%%%%%% VERBOSE:		Displays run messages
%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: Jul 11 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AtomsToStats(in_dir, ff_id, out_dir, L, ipr_factor, angle_th, grain_th, RMN1, RMANG, rmang_th, PFIG, PFIG_N, VERBOSE)
init_time = tic;

%This code makes use of random integers to prevent spatial biasing
rng('shuffle');
rseed = rng; rseed = rseed.Seed;
if VERBOSE fprintf('Initialized with random seed %d.\n', rseed); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for fid = ff_id
	%Initialize time and filenames
	start_time = tic;
	fname = sprintf('%sAtomicPos_%d.mat', in_dir, fid);
	if ~exist(fname, 'file')
		break
	end
	
	%Import the atomic positions
	load(fname);
	if VERBOSE fprintf('Loaded file ID %d.\n', fid); end
	
	%Compute the interatomic distance from the domain size and the number of atoms
	%This assumes that all atoms have the same hexagonal unit cell and equally share the SQUARE domain
	%The method should work for rectangular domains, but this needs to be changed for the domain's area
	interatomic_distance = sqrt(2*(L^2/length(atoms_x))/sqrt(3));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Pad the atomic positions with ghosts
	sub_time = tic;
	[atoms_x, atoms_y, ghost_ind, n_true, n_fake] = MakeGhosts(atoms_x, atoms_y, L, interatomic_distance);
	
	%Obtain the Voronoi tesselation to obtain the cells and their polygon's vertices
	[voronoi_v, voronoi_c] = voronoin([atoms_x', atoms_y']);
	
	%Now compute the Voronoi cells' properties
	atoms_angle = zeros(1, n_fake);
	atoms_area = atoms_angle; atoms_perim = atoms_angle;
	for cur_atom = 1:n_fake
		%Extract the vertices of the current cell's polygon
		poly_x = voronoi_v(voronoi_c{cur_atom}', 1); poly_y = voronoi_v(voronoi_c{cur_atom}', 2);
		
		%Compute the area and perimeter
		spoly_x = [poly_x(end); poly_x(1:end-1)]; spoly_y = [poly_y(end); poly_y(1:end-1)];
		atoms_area(cur_atom) = 0.5*abs(sum(poly_x .* spoly_y - poly_y .* spoly_x));
		atoms_perim(cur_atom) = sum(sqrt((poly_x-spoly_x).^2 + (poly_y-spoly_y).^2));
		
		%Compute the average angle between the atom and the vertices
		atoms_angle(cur_atom) = CircMean(mod(atan2(atoms_y(cur_atom)-poly_y, atoms_x(cur_atom)-poly_x), pi/3), pi/3);
	end
	
	%Compute the measure of niceness to help separate grains
	bad_atoms = ipr_factor * abs(IPR(atoms_area, atoms_perim) - (pi/6)/tan(pi/6));
	
	%Replace NaNs with a large number because Matlab is dumb
	bad_atoms(isnan(bad_atoms)) = max(bad_atoms);
	bad_atoms = (bad_atoms >= 1.0);
	
	if VERBOSE fprintf('Computed Voronoi properties in %.1fs.\n', toc(sub_time)); end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Obtain the adjacency matrix of the atoms
	sub_time = tic;
	adj_matrix = PointAdjMatrix(atoms_x, atoms_y);
	if VERBOSE fprintf('Computed adjacency matrix in %.1fs.\n', toc(sub_time)); end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Find the grains using a flood filling algorithm 
	sub_time = tic;
	[grain_n, atoms_grain_id, g_angle] = FindGrains(n_true, n_fake, atoms_angle, ghost_ind, bad_atoms, ...
		adj_matrix, grain_th, angle_th);
	if VERBOSE fprintf('Done flooding grains. Now filling the gaps.\nFound '); end
	
	%We now fill in the gaps between grains
	temp_grain_id = atoms_grain_id;
	while true
		%We must convert true and ghost atoms since ghosts are also used later on
		fetch_id = find(atoms_grain_id <= 0);
		if isempty(fetch_id) break; end
		if VERBOSE fprintf('%d, ', length(fetch_id)); end
		
		for check_i = 1:length(fetch_id)
			%Pick the neighbors's ids then pick their grain ids
			neighbors = adj_matrix{fetch_id(check_i)};
			
			%Remove neighbors that are not yet assigned
			neighbors = neighbors(atoms_grain_id(neighbors) > 0);
			
			%If no neighbor is close to a grain, continue
			if isempty(neighbors)
				continue;
			end
			
			%Pick the mode of the valid grain ids. Mode only returns the first id it finds, so the method might
			%be biased towards lower id grains and hence larger grains. This doesn't seem to affect the results.
			assign_grain = -mode(-atoms_grain_id(neighbors));
			
			%Assign the mode to the current atom
			temp_grain_id(fetch_id(check_i)) = assign_grain;
		end
		
		atoms_grain_id = temp_grain_id;
	end
	if VERBOSE fprintf('unassigned atoms.\nDetected %d grains in %.1fs.\n', grain_n, toc(sub_time)); end
	
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%We now find grain boundaries and measure grain properties for the first time
	sub_time = tic;
	[g_area, g_perim, g_adjmat, int_length, mis_angle] = MeasureGrains(g_angle, ...
				atoms_grain_id, ghost_ind, n_true, grain_n, adj_matrix, atoms_area(1:n_true), interatomic_distance);
	
	if VERBOSE
		fprintf('Computed grain properties in %.1fs.\n', toc(sub_time));
		fprintf('Averages... Area: %.1f, Perimeter: %.1f.\n', mean(g_area), mean(g_perim));
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Post process the grain distribution now that we have their properties to remove spurrious grains
	if RMN1 || RMANG
%		%The loop can be implemented to remove all grains with a single neighbor or small misorientations
%		%RECURSIVELY guaranteeing that these filters are sharp even for very fine distributions
% 		while true
			sub_time = tic;
			
			%Identify grains to be removed because they are inside others
			if RMN1
				delete_id = find(sum(g_adjmat, 2) == 1);
				parent_id = delete_id;
				for traverse_id = 1:length(delete_id)
					parent_id(traverse_id) = find(g_adjmat(delete_id(traverse_id), :));
				end
			else
				delete_id = [];
				parent_id = [];
			end
			
			%Identify grains to be removed because their angle is too similar to a neighbor
			if RMANG
				check_ij = find(triu(g_adjmat));
				[g_i, g_j] = ind2sub(size(g_adjmat), check_ij(mis_angle(check_ij) < rmang_th));
				delete_gi = (g_area(g_i) <= g_area(g_j));
				delete_id = [delete_id; delete_gi.*g_i + (1-delete_gi).*g_j];
				parent_id = [parent_id; delete_gi.*g_j + (1-delete_gi).*g_i];
			end
			
			%Filter out repetitions in case a grain is excluded for multiple reasons
			[delete_id, order] = unique(delete_id);
			parent_id = parent_id(order);
			
			if ~isempty(delete_id)
				%Loop over grains to be deleted
				for tt = 1:length(delete_id)
					%We now replace the atomic IDs of the current grain by the parent's
					atoms_grain_id(atoms_grain_id == delete_id(tt)) = parent_id(tt);
					parent_id(parent_id == delete_id(tt)) = parent_id(tt);
				end
				
				%We must restore the linear 1:n structure of grain IDs
				old_ids = unique(atoms_grain_id);
				old_grain_n = grain_n; grain_n = length(old_ids);
				temp_grain_id = atoms_grain_id;
				for tt = 1:grain_n
					atoms_grain_id(temp_grain_id == old_ids(tt)) = tt;
				end
				
				%We must also rework the angle matrix
				g_angle = g_angle(old_ids);
				
				%Finally recompute the grain properties after the removal
				[g_area, g_perim, g_adjmat, int_length, mis_angle] = MeasureGrains(g_angle, ...
					atoms_grain_id, ghost_ind, n_true, grain_n, adj_matrix, atoms_area(1:n_true), interatomic_distance);
				
				if VERBOSE
					fprintf('Recomputed grain properties after post-processing in %.1fs.\n', toc(sub_time));
					fprintf('Deleted %d grains.\n', old_grain_n - grain_n);
				end
% 			else
% 				break
% 			end
		end
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Now extract the convex hull measure of regularity
	%These correspond to the regions defined by the atomic positions, not the full grains. The convex hull
	%ratio will have to be corrected with the perimeter times the boundary thickness.
	sub_time = tic;
	[g_hull_area, g_hull_perim] = MeasureConvexHull(atoms_grain_id, atoms_x, atoms_y, n_true, grain_n);
	if VERBOSE
		fprintf('Computed convex hulls in %.1fs.\n', toc(sub_time));
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Print figures
	if PFIG
		PrintAng(sprintf('%sAng_%d_R.png', out_dir, fid), n_true, atoms_angle, atoms_x, atoms_y, PFIG_N, L);
		PrintAng(sprintf('%sAng_%d.png', out_dir, fid), n_true, g_angle(atoms_grain_id), atoms_x, atoms_y, PFIG_N, L);
	end
	
	%Save the processed variables into the output file
	fprintf('Converted file %d in %.1fs.\n', fid, toc(start_time));
	save(sprintf('%sData_%d.mat', out_dir, fid), 'L', 'n_true', 'grain_n', ...
		'g_angle', 'g_area', 'g_perim', 'g_adjmat', 'int_length', 'mis_angle', 'g_hull_area', 'g_hull_perim');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Print the elapsed time and exit
total_time = toc(init_time);
fprintf('Time elapsed in AtomsToStats: %dh %dm %ds\n', floor(total_time/3600), ...
	floor(mod(total_time, 3600)/60), ceil(mod(total_time, 60)))
end
