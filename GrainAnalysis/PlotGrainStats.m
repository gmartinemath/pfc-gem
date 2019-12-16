%Plots grain statistics
%Gabriel Martine
%Last updated: May 16 2019

function PlotGrainStats(stats, area_to_diameter, min_coord, max_coord, color, thickness, mark)
	%Plot the diameter if needed
	if area_to_diameter
		[x, n, error] = HistMaker_FD(sqrt(stats.area)/mean(sqrt(stats.area)), true);
	else
		[x, n, error] = HistMaker_FD(stats.area/mean(stats.area), true);
	end
	figure(1); errorbar(x, n, error, color, 'linewidth', thickness);
	
	%Normalize the area and perimeter by their average
	stats.area = stats.area / mean(stats.area);
	stats.perim = stats.perim / mean(stats.perim);
		
	%Simple statistics and CHR
	%figure(2); [x, n] = HistMaker_FD(stats.perim, true); plot(x, n, color, 'linewidth', thickness);
	%figure(3); [x, n] = HistMaker_FD(stats.ipr, true); plot(x, n, color, 'linewidth', thickness);
	%figure(4); [x, n] = HistMaker_FD(stats.chr, true); plot(x, n, color, 'linewidth', thickness);

	figure(2); [x, n, error] = HistMaker_FD(stats.perim, true); errorbar(x, n, error, color, 'linewidth', thickness);
	figure(3); [x, n, error] = HistMaker_FD(stats.ipr, true); errorbar(x, n, error, color, 'linewidth', thickness);
	figure(4); [x, n, error] = HistMaker_FD(stats.chr, true); errorbar(x, n, error, color, 'linewidth', thickness);
	
	%Statistics wtr coordination number
	coord = [];
	area_vs_coord = [];
	error_area_vs_coord = [];
	avg_nc_vs_coord = [];
	error_avg_nc_vs_coord = [];
	for cc = min_coord:max_coord
		if any(stats.coord == cc)
			coord(end+1) = sum(stats.coord == cc);
			area_vs_coord(end+1) = mean(stats.area(stats.coord == cc));
			error_area_vs_coord(end+1) = std(stats.area(stats.coord == cc))/sqrt(sum(stats.coord == cc));
			avg_nc_vs_coord(end+1) = mean(stats.avg_neighbor_coord(stats.coord == cc));
			error_avg_nc_vs_coord(end+1) = std(stats.avg_neighbor_coord(stats.coord == cc))/sqrt(sum(stats.coord == cc));
		else
			coord(end+1) = 0;
			area_vs_coord(end+1) = 0; error_area_vs_coord(end+1) = 0;
			avg_nc_vs_coord(end+1) = 0; error_avg_nc_vs_coord(end+1) = 0;
		end
	end
	
	%Histogram, so error is just sqrt(coord)
	figure(5); errorbar(min_coord:max_coord, coord/trapz(coord), sqrt(coord)/trapz(coord), color, 'linewidth', thickness);%, 'marker', mark, 'markersize', 6, 'markerfacecolor', color);
	
	%Value is an average so error is the standard error sigma_n/sqrt(N_n)
	figure(6); errorbar(min_coord:max_coord, area_vs_coord, error_area_vs_coord, color, 'linewidth', thickness);%, 'marker', mark, 'markersize', 6, 'markerfacecolor', color);
	
	%Error propagation using the previous errors
	error_combo = coord.*area_vs_coord .* sqrt((sqrt(coord)./coord).^2 + (error_area_vs_coord./area_vs_coord).^2);
	figure(7); errorbar(min_coord:max_coord, coord.*area_vs_coord/trapz(coord.*area_vs_coord), error_combo/trapz(coord.*area_vs_coord), color, ...
		'linewidth', thickness);%, 'marker', mark, 'markersize', 6, 'markerfacecolor', color);
	
	%Value is an average so error is the standard error sigma_n/sqrt(N_n)
	figure(8); errorbar(min_coord:max_coord, avg_nc_vs_coord, error_avg_nc_vs_coord, color, 'linewidth', thickness);%, 'marker', mark, 'markersize', 6, 'markerfacecolor', color);
end
