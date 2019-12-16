%This automatically fills in the histogram parameters by choosing the bin number according to the Friedman-Diaconis rule
%The left and right bounds are chosen using the extrema of the distribution
%Gabriel Martine
%Last updated: Oct 10 2018

function [grid, histg, stdg] = HistMaker_FD(data, norm_freq)
	bin_count = (max(data)-min(data))/(2*iqr(data)/length(data)^(1/3));	
	[grid, histg, stdg] = HistMaker(data, min(data), max(data), bin_count, norm_freq);
end
