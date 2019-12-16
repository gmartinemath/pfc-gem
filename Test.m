%PFC parameters
n_atoms = 256;			%The number of atoms in each dimension
ppa = 8;				%The number of pixels per atoms, this should be at least 8
m = 0.07; beta = 0.025;	%PFC parameters
atom_level = -0.035;	%Extraction level for atoms, depends on (m, beta)

%Run the PFC evolution
pfc_run_output_dir = '~/Downloads/PFC_Test/PFC_Run/G1/';
mkdir(pfc_run_output_dir);
PFC_Evolution('', pfc_run_output_dir, 'NOISE', ...	%Use a random initial condition
	n_atoms, ppa, ...								%Domain size
	1, 1000, 10, 10, ...							%Runs from time 1 to 1000, saves 10 exponentially spaced snapshots
	0.05, ...										%Scale of random noise in the initial condition
	atom_level, m, beta, ...						%Extraction level and PFC parameters
	600, 2*beta, ...								%Scheme parameters: tau = 1000 and C is always taken as 2*beta
	false);											%Verbose output


%Atom extraction is done during the simulation to minimize the memory footprint


%Measurement parameters
load('~/Downloads/PFC_Test/PFC_Run/G1/MetaData.mat', 'L');	%This fetches the physical domain size
gamma = 0.001;												%Isoperimetric ratio closeness check
theta = 2.5*pi/180;											%Threshold misorientation
alpha = 40;													%Threshold number of atoms per grain

%Run the grain and measurement extraction
ats_output_dir = '~/Downloads/PFC_Test/ATS_Run/G1/';
mkdir(ats_output_dir);

AtomsToStats(pfc_run_output_dir, [9, 10], ats_output_dir, ...	%Specify to extract the last two atom distributions
	L, 1/gamma, theta, alpha, ...								%Processing parameters
	true, ...													%Post-processing: Remove grains with only one neighbor
	true, theta, ...											%Post-processing: Remove grains with misorientation below theta
	true, n_atoms, ...											%Save colored pictures with resolution n_atoms
	false);														%Verbose output


%Start the statistics plotting
OverallAnalysis
