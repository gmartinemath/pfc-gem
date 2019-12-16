%This function pads the given atomic distribution with ghosts to allow periodic computations
%Gabriel Martine
%Last updated: July 17 2017

function [atoms_x, atoms_y, ghost_ind, n_true_atoms, n_fake_atoms] = MakeGhosts(atoms_x, atoms_y, L, pad_dist)
	
	%We must pad the atomic positions to compute the periodic Voronoi and Delaunay tesselations. However, the
	%ghost atoms must still "point" to the original atoms for the adjacency matrix and graph calculations
	n_true_atoms = length(atoms_x);
	
	%Prepare an index map which is the identity over the true atoms
	ghost_ind = (1:n_true_atoms);
	
	%Find all atoms that are close to the horizontal boundaries and replicate them to the other side
	roll_low = find(atoms_x < 2*pad_dist); roll_high = find(atoms_x >= L-2*pad_dist);
	ghost_ind = [ghost_ind, ghost_ind(roll_low), ghost_ind(roll_high)];
	atoms_x = [atoms_x, atoms_x(roll_low), atoms_x(roll_high)];
	atoms_y = [atoms_y, atoms_y(roll_low), atoms_y(roll_high)];
	shifts_x = [zeros(1, n_true_atoms), -L*ones(1, length(roll_low)), L*ones(1, length(roll_high))];
	shifts_y = 0*shifts_x;
	
	%Repeat the above for the vertical boundaries, including the new ghost atoms
	roll_low = find(atoms_y < 2*pad_dist); roll_high = find(atoms_y >= L-2*pad_dist);
	ghost_ind = [ghost_ind, ghost_ind(roll_low), ghost_ind(roll_high)];
	atoms_x = [atoms_x, atoms_x(roll_low), atoms_x(roll_high)];
	atoms_y = [atoms_y, atoms_y(roll_low), atoms_y(roll_high)];
	shifts_x = [shifts_x, shifts_x(roll_low), shifts_x(roll_high)];
	shifts_y = [shifts_y, -L*ones(1, length(roll_low)), L*ones(1, length(roll_high))];
	
	atoms_x = atoms_x - shifts_x; atoms_y = atoms_y - shifts_y;
	
	%Now all valid atoms are within [0, L)x[0, L) and consist in the n_true_atoms first entries in the arrays
	n_fake_atoms = length(atoms_x);
end
