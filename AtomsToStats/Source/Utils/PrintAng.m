%This function prints the angle image of atoms at (x,y) at a given angle 
%Gabriel Martine
%Last updated: May 8 2017

function PrintAng(fname, n_ang, angle, x, y, N, L)
	%Create arrays to accumulate angles
	buff_grains = -ones(N);
	
	%Accumulate angles onto the grid, ignoring repetitions. This does not matter since inside grains, the repetitions
	%will be close while repetitions at the boundary are essentially randomized anyways
	bx = floor(mod(x, L)*N/L)+1; by = floor(mod(y, L)*N/L)+1;
	for i=1:n_ang
		buff_grains(by(i), bx(i)) = angle(i);
	end
	
	%Replace undefined pixels with NaN
	buff_grains(buff_grains == -1) = NaN;
	
	%Save the new grain image
	DumpImage(buff_grains, fname, 'hsv', 0.0, pi/3, 256);
end

