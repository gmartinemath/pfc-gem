%Usual PFC parameters
beta = 0.025; m = 0.07; phase_scale = 0.05; atom_level = -0.035; n_atoms = 256; ppa = 32; big_tau = 1000;

for GG = [1,2]
%Create the initial condition
t_initial = 0; t_final = 40; snapshots = 0;
PFC_Evolution('', sprintf('~/MemGEM/InitialConditions/PFC_Standard_256x32/G%d/',GG), ...
	'NOISE', n_atoms, ppa, t_initial, t_final, snapshots, snapshots, phase_scale, atom_level, m, beta, big_tau, 2*beta);

%Evolve from 40 to 40000
t_initial = 40; t_final = 20047; snapshots = 45;
PFC_Evolution(sprintf('~/MemGEM/InitialConditions/PFC_Standard_256x32/G%d/FinalFull.mat',GG)', sprintf('~/MemGEM/PFC_Standard_256x32/G%d/',GG), ...
	'PHASE', n_atoms, ppa, t_initial, t_final, snapshots, snapshots, phase_scale, atom_level, m, beta, big_tau, 2*beta);
end
