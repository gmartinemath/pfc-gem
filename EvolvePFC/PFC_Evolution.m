%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% PFC EVOLUTION
%%%%%%%%%%%%%% This code will load an initial condition in the form of a matlab array
%%%%%%%%%%%%%% with values corresponding to a lattice angle. This is used to create an
%%%%%%%%%%%%%% initial PFC phase distribution which is evolved in time. During the
%%%%%%%%%%%%%% evolution, the atomic positions are saved periodically for further analysis.
%%%%%%%%%%%%%% The code will also save arrays containing the PFC energy of the phase. 
%%%%%%%%%%%%%% 
%%%%%%%%%%%%%% CLI INPUT VARIABLES
%%%%%%%%%%%%%% in_file:		The filename of the input phase or file containing the random seed
%%%%%%%%%%%%%% out_dir:		The directory where the output (.mat) should be saved
%%%%%%%%%%%%%% init_cond:	The type of initial conditions: 'NOISE' for a noisy initial
%%%%%%%%%%%%%%				condition or 'PHASE' if the initial condition is a PFC phase
%%%%%%%%%%%%%% n_atoms: 	The (approximate) number of atoms in each dimension
%%%%%%%%%%%%%% ppa:			The (approximate) number of pixels per atom in each dimension
%%%%%%%%%%%%%% begin_ind:	Given the same tau and C, the time step of the initial condition
%%%%%%%%%%%%%% end_ind:		The time step of the final array that is computed
%%%%%%%%%%%%%% snap_energy:	The number of energy datapoints to be computed
%%%%%%%%%%%%%% snap_atoms:	The number of atomic position snapshots to be saved for further analysis
%%%%%%%%%%%%%% phase_scale: The scaling of the initial condition
%%%%%%%%%%%%%% height:		The level at which to extract atoms
%%%%%%%%%%%%%% m:			The average mass
%%%%%%%%%%%%%% beta:		The inverse temperature
%%%%%%%%%%%%%% tau:			The physical time step size
%%%%%%%%%%%%%% C:			The regularization constant
%%%%%%%%%%%%%% VERBOSE:		Displays run messages
%%%%%%%%%%%%%% 
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: Dec 16 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PFC_Evolution(in_file, out_dir, init_cond, n_atoms, ppa, begin_ind, end_ind, ...
	snap_energy, snap_atoms, phase_scale, height, m, beta, tau, C, VERBOSE)

init_time = tic;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Physical domain size, the number of grid points and the physical grid
[Nx, L, gridx, hx] = PrepareGrid_PFC(n_atoms, ppa);
Ny = Nx; gridy = gridx; hy = hx; %Legacy! Didn't test rectangular domain in ages so they may not work

%Time indices at which the energy and atomic positions are saved
save_energy = SaveTimes(begin_ind, end_ind, snap_energy);
save_atoms = SaveTimes(begin_ind, end_ind, snap_atoms);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rseed = 0;
switch init_cond
	case 'PHASE'
		load(in_file);
		if ~(exist('cur_phase', 'var'))
			fprintf('Cannot load initial condition.\n');
			return;
		end
	case 'NOISE'
		%The initial condition is random noise
		if exist(in_file, 'file')
			load(in_file, 'rseed');
			rng(rseed);
		else
			rng('shuffle');
			rseed = rng; rseed = rseed.Seed;
			if VERBOSE fprintf('\nInitialized with random seed %d.\n', rseed); end
		end
				
		cur_phase = phase_scale*randn(Ny, Nx);
		cur_phase = cur_phase + m - mean(cur_phase(:));
		save(sprintf('%sInitialFull.mat', out_dir), 'cur_phase','-v7.3'); 
	otherwise
		fprintf('Invalid initial condition.\n');
		return;
end

%Save the input values and phase for further processing
save(sprintf('%sMetaData.mat', out_dir), 'n_atoms', 'ppa', 'm', 'beta', 'tau', 'C', 'L', 'phase_scale', 'rseed');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating the Fourier differentiation operators
[X, Y] = meshgrid(0:Nx-1, 0:Ny-1);
F = 2.0/hx^2 * (cos(2.0*pi*X/Nx) - 1.0) + 2.0/hy^2 * (cos(2.0*pi*Y/Ny) - 1.0);
G_PFC = (1.0 - tau * F .* ((F+1.0).^2 + C - beta)).^(-1);

%Cleaning up
clearvars X Y;

%Initializing working arrays
fft_phase = fft2(cur_phase);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up arrays to collect the energy
energy_PFC = 10e10*ones(1, length(save_energy));
times_energy = energy_PFC;
times_atoms = 10e10*ones(1, length(save_atoms));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin the iteration
ind_t = begin_ind; ind_array = 1; ind_energy = 1;
cur_time = tic;
if VERBOSE fprintf('Beginning iteration.\n'); end
while (ind_t < end_ind)
	%Update the PFC phase with the original equation
	fft_phase = G_PFC .* (fft_phase + tau * F .* (fft2(cur_phase.^3) - C*fft_phase));
	
	cur_phase = real(ifft2(fft_phase));
	fft_phase = fft2(cur_phase);
	
	ind_t = ind_t + 1;
	
	
	%Save energy once in a while
	if save_energy(ind_energy) == ind_t
		%PFC Energy according to the PFC functional
		energy_PFC(ind_energy) = mean2(0.5*cur_phase .* real(ifft2((F+1.0).^2 .* fft_phase)) ...
			+ 0.25*(cur_phase.^2 - beta).^2);
		times_energy(ind_energy) = ind_t;
		
		%Print the time taken in the last energy_skip steps
		if VERBOSE fprintf('Saved energy at step %d (spent %.1fs since last print).\n', ind_t, toc(cur_time)); end
		cur_time = tic;
		
		%Save the current energy arrays just in case the simulation crashes
		save(sprintf('%sEnergy_PFC.mat', out_dir), 'energy_PFC');
		save(sprintf('%sTimes_Energies.mat', out_dir), 'times_energy');
		ind_energy = ind_energy + 1;
	end
	
	%Save atomic positions every once in a while
	if save_atoms(ind_array) == ind_t
		[atoms_x, atoms_y] = FindAtoms(cur_phase, height, ppa, L, L);
		save(sprintf('%sAtomicPos_%d.mat', out_dir, ind_array), 'atoms_x', 'atoms_y');
		
		%UNCOMMENT TO ALSO SAVE THE PHASE ITSELF
		%save(sprintf('%sAtomicPhase_%d.mat', out_dir, ind_array), 'cur_phase');
		
		times_atoms(ind_array) = ind_t;
		
		%Save the current time array just in case the simulation crashes
		save(sprintf('%sTimes_Atoms.mat', out_dir), 'times_atoms');
		ind_array = ind_array + 1;
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save the full-size array in case we need to continue the calculations
save(sprintf('%sFinalFull.mat', out_dir), 'cur_phase','-v7.3');

%Save the final time and energy arrays
save(sprintf('%sEnergy_PFC.mat', out_dir), 'energy_PFC');
save(sprintf('%sTimes_Energies.mat', out_dir), 'times_energy');
save(sprintf('%sTimes_Atoms.mat', out_dir), 'times_atoms')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Print the elapsed time and exit
total_time = toc(init_time);
fprintf('Time elapsed in PFC_Evolution: %dh %dm %ds\n', floor(total_time/3600), ...
	floor(mod(total_time, 3600)/60), ceil(mod(total_time, 60)));
end
