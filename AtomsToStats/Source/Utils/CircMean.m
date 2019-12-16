%Computes the angular mean of data in the periodic interval [0, L). This is the only way to have a consistent periodic
%mean. If dat is sufficiently clustered, this mean is consistent with the usual mean on R
%Gabriel Martine
%Last updated: Jan 3 2017

function avg = CircMean(dat, L)
	rfac = 2.0*pi/L;
	avg = atan2(mean(sin(dat*rfac)), mean(cos(dat*rfac)))/rfac;
	avg = mod(avg, L);
end
