%Creates a grid and the histogram values of data in between the values x_l and x_r, discarding anything out of the bounds
%The usual histogram function can be recovered with [x, n] = HistMaker(dat, min(dat), max(dat), n_bin, false);
%Thus, hist(dat, n_bin); plot(x, n); will coincide
%Gabriel Martine
%Last updated: Oct 10 2018

function [grid, histg, stdg] = HistMaker(data, x_l, x_r, n_bin, norm_freq)
	n_bin = round(n_bin);
	
	%Build a grid from x_l to x_r with n_bin points
	grid = linspace(x_l, x_r, n_bin+1); grid(end) = [];
	bin_width = grid(2) - grid(1);
	
	%Bin the data onto the grid using a simple [,) rule
	histg = zeros(1, n_bin);
	stdg = zeros(1, n_bin);
	for i = 1:n_bin
		binned_id = find((data >= x_l + bin_width*(i-1)) & (data < x_l + bin_width*i));
		histg(i) = length(binned_id);
		
		%IS THIS CORRECT?
		stdg(i) = sqrt(length(binned_id));	%Looks good - assumes Poisson process...
	end
	
	%Shift the x values so they lie in the middle of the bins
	grid = grid + bin_width/2;
	
	%Make the histogram into a probably density distribution
	if norm_freq
		hist_area = trapz(grid, histg);
		histg = histg / hist_area;
		stdg = stdg / hist_area;
	end
end
