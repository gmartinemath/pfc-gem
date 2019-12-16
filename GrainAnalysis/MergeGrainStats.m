%Utility to merge two distribution statistics
%Gabriel Martine
%Last updated: Apr 10 2019

function merged_stats = MergeGrainStats(stats1, stats2)
	if ~isstruct(stats1)
		merged_stats = stats2;
		return
	end
	
	merged_stats.area = [stats1.area; stats2.area];
	merged_stats.perim = [stats1.perim; stats2.perim];
	merged_stats.ipr = [stats1.ipr; stats2.ipr];
	merged_stats.chr = [stats1.chr; stats2.chr];
	merged_stats.coord = [stats1.coord; stats2.coord];	
	merged_stats.avg_neighbor_coord = [stats1.avg_neighbor_coord; stats2.avg_neighbor_coord];
end
