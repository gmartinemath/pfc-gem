%Computes the angle distance between two angles in [0, L]. Note that the maximum output of ad is L/2
%Gabriel Martine
%Last updated: Apr 7 2017

function ad = AngleDist(d, L)
	ad = abs(mod(d + L/2, L) - L/2);
end
