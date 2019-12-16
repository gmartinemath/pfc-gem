%Files
dirs = {'PFC_PhaseScale001', 'PFC_PhaseScale010', 'PFC_V1_1024x8', 'PFC_V2_1024x8', 'PFC_V4_1024x8', 'PFC_V5_1024x8'};
run_id = 1:3;
ff_id = [1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50];

%Parameters
ang = 2.5*pi/180;
arr = 40;
gamma = 1000;
for cur_dir = 1:length(dirs)
	load(sprintf('~/MemGEM/%s/G%d/', dirs{cur_dir}, run_id(cur_run)), 'L', 'n_atoms');
	for cur_run = 1:length(run_id)
		mkdir(sprintf('~/MemGEM/ExtractedConvexHull/%s/G%d/', dirs{cur_dir}, run_id(cur_run)));
		AtomsToStats(sprintf('~/MemGEM/%s/G%d/', dirs{cur_dir}, run_id(cur_run)), ff_id, ...
			sprintf('~/MemGEM/ExtractedConvexHull/%s/G%d/', dirs{cur_dir}, run_id(cur_run)), ...
			L, gamma, ang, arr, true, true, ang, true, n_atoms, false);
	end
end
