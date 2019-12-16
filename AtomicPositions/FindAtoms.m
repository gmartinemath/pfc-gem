%Find atoms in a periodic phase at the given height level.
%ppa is the approximate number of pixels per atoms, used to extend the image and allow atom detection near boundaries
%Lx and Ly are the boundary lengths but note that the subsequent grain analysis assumes L = Lx = Ly for simplicity
%Note that the domain must be periodic so Lx and Ly do NOT correspond to the very last entry in phase. We have that
%phase(1,1) is at (0,0) and phase(end,end) is at (Lx-hx,Ly-hy) since (Lx,Ly) is identified to (0,0)
%Gabriel Martine
%Last updated: Jul 11 2018

function [atoms_x, atoms_y] = FindAtoms(phase, level, ppa, Lx, Ly)
	%Defining the grid so the results can be reported in the domain [0,Lx) x [0,Ly)
	Nx = size(phase, 2); gridx = linspace(0, Lx, Nx+1); gridx(end) = [];
	Ny = size(phase, 1); gridy = linspace(0, Ly, Ny+1); gridy(end) = [];
	
	%The actual maximum number of atoms in each direction. 16 atoms is a good value
	atoms_p_img = 16;
	ppa = ceil(ppa); %We need an integer ppa so take the ceiling
	m_pix = atoms_p_img*ppa;
	
	%A zero array used to accumulate atoms every loop instead of growing arrays
	%In general, we have about 1.45*n_atoms^2 but we should pick a bit more
	expected_array = zeros(1, round(1.6*atoms_p_img^2));
	
	%Pad the phase and the grids to allow for periodic detection of contours. The grids
	%are extended as non-periodic reals to allow cropping the atoms
	phase_pad = wextend(2, 'ppd', phase, ppa);
	gridx_pad = wextend(1, 'sp1', gridx, ppa);
	gridy_pad = wextend(1, 'sp1', gridy, ppa);
	
	%Find out how many times the image must be split and compute how much would spill onto a new subimage
	Wy = floor(Ny/m_pix); Wx = floor(Nx/m_pix);
	remainder_y = Ny - Wy*m_pix; remainder_x = Nx - Wx*m_pix;
	
	%Loop over the padded subimages taking care to change the endpoint of the leftover region
	atoms_x = []; atoms_y = [];
	for iy = 0:Wy
		ind_y1 = 1 + iy*m_pix;
		ind_y2 = (iy+1)*m_pix + 2*ppa;
		if iy == Wy
			if remainder_y == 0 break; end
			ind_y2 = Ny + 2*ppa;
		end
		
		for ix = 0:Wx
			ind_x1 = 1 + ix*m_pix;
			ind_x2 = (ix+1)*m_pix + 2*ppa;
			if ix == Wx
				if remainder_x == 0 break; end
				ind_x2 = Nx + 2*ppa;
			end
			
			%Compute the atomic positions on the subimage
			[p_x, p_y] = GetContourCentroids(phase_pad(ind_y1:ind_y2, ind_x1:ind_x2), level, ...
				gridx_pad(ind_x1:ind_x2), gridy_pad(ind_y1:ind_y2), expected_array);
			
			%Add the partial lists while removing false atoms, those outside of the padded subimage
			out_bounds = (p_y < gridy_pad(ind_y1+ppa) | p_y >= gridy_pad(ind_y2-ppa+1) | ...
				p_x < gridx_pad(ind_x1+ppa) | p_x >= gridx_pad(ind_x2-ppa+1));
			
			atoms_x = [atoms_x, p_x(~out_bounds)];
			atoms_y = [atoms_y, p_y(~out_bounds)];
		end
	end
end
