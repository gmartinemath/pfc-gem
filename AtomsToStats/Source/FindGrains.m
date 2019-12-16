%This will assign an id to each atom depending on which grain they are a part of
%A value of 0 means the atom is unassigned. A value of -1 means the atom is rejected
%as a starting point because it is not a nice atom. A value of -2 means the atom is rejected
%as a starting point because it was part of a too small grain.
%Gabriel Martine
%Last updated: May 11 2017

function [grain_n, atoms_grain_id, grain_angle] = FindGrains(n_true, n_fake, angles, ghost_ind, ...
	bad_atoms, adj_matrix, grain_th, angle_th)
	
	%Initialize the grain ID vector with -1 where atoms are rejected because they are not nice enough
	atoms_grain_id = -bad_atoms;
	
	%Initialize grain number and angles
	grain_n = 0;
	grain_angle = [];
	
	%Get the atoms that have a ghost
	atom_has_ghost = false(1, n_fake);
	atom_has_ghost(ghost_ind(n_true+1:end)) = true;
	
	%Loop over the unassigned nice atoms and flood fill their grain
	while true
		%To prevent directional bias, I absolutely want some sort of randomization of position...
		initial_ind = find(atoms_grain_id(1:n_true) == 0);
		
		%If there are no valid starting atoms left, break the loop
		if isempty(initial_ind)
			break;
		end
		
		%Otherwise, pick a random id
		initial_ind = initial_ind(randi(length(initial_ind)));
		
		%Initialize "dangling" atoms at the active edge of the growing flood 
		dangling = initial_ind;
		temp_dangling = [];
		
		%Initialize mask arrays to capture the current grain and its ghosts
		grain_msk = false(1, n_fake);
		ghost_msk = grain_msk;
		grain_msk(initial_ind) = true;
		
		%Create a copy of atoms_grain_id since it only needs it to be false where good atoms can still be filled
		valid_atoms = logical(~atoms_grain_id);
		
		cur_angle = angles(initial_ind);
		%While the boundary of the flood is still active
		while ~isempty(dangling)
			for check_i = dangling
				for check_j = adj_matrix{check_i}
					%For each neighbor of the current dangling atom, check their niceness
					if valid_atoms(check_j) && ~grain_msk(check_j) && ...
							(AngleDist(angles(check_j) - cur_angle, pi/3) < angle_th)
						%Get the true atom's index
						original_j = ghost_ind(check_j);
						
						%Add the true atom to the grain mask and the dangling list
						grain_msk(original_j) = true;
						temp_dangling = [temp_dangling, original_j];
						
						%Deactivate the parent atom's ghosts
						if atom_has_ghost(check_i)
							check_ghost = find((ghost_ind == check_i));
							check_ghost = check_ghost(2:end);
							valid_atoms(check_ghost) = false;
							ghost_msk(check_ghost) = true;
						end
						
						%Deactivate the new atom's ghosts
						if atom_has_ghost(original_j)
							check_ghost = find((ghost_ind == original_j));
							check_ghost = check_ghost(2:end);
							valid_atoms(check_ghost) = false;
							ghost_msk(check_ghost) = true;
						end
					end
				end
			end
			
			%Compute the mean angle of the currently grown grain 
			cur_angle = CircMean(angles(grain_msk), pi/3);
			
			%Replace the ancient dangling atoms by the new boundary layer
			dangling = temp_dangling;
			temp_dangling = [];
		end
		
		
		%Check if the grain is large enough
		if sum(grain_msk) < grain_th
			%If not, mark the detected atoms as bad
			atoms_grain_id(grain_msk) = -2;
			atoms_grain_id(ghost_msk) = -2;
		else
			%If it is large enough, save the grain
			grain_n = grain_n + 1;
			
			atoms_grain_id(grain_msk) = grain_n;
			grain_angle = [grain_angle; cur_angle];
			
			%This is used to make the grain assignment and boundary determination accurate at the domain's boundary
			atoms_grain_id(ghost_msk) = grain_n;
		end
	end
end
