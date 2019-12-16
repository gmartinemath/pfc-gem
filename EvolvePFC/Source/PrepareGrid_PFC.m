%Outputs the PFC domain
%n_atoms is the number of atoms in the lattice direction aligned with an axis
%pix_p_atoms or ppa is the number of pixels per atoms in a direction aligned with an axis
%These values must be integers and should be powers of 2
%Gabriel Martine
%Last updated: May 4 2017

function [N, L, grid, h] = PrepareGrid_PFC(n_atoms, pix_p_atoms)
	L = 4.0*pi/sqrt(3) * n_atoms; N = n_atoms * pix_p_atoms;
	grid = linspace(0,L,N+1); grid(end) = []; h = grid(2)-grid(1);
end
