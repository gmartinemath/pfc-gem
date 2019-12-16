%Find contour centroids and associated properties in a given image
%Gabriel Martine
%Last updated: May 4 2017

function [atoms_x, atoms_y] = GetContourCentroids(img, level, gridx, gridy, exp_n)
	%Compute the level-contour of img using the given 2d grids. The output cc is a 2x(S+T) matrix where
	%S is the number of connected curves and T is the total number of points identified in those curves.
	%For every connected curve, the matrix has an entry (level, t) where t is the number of points on
	%the current curve. The following t entries are of the form (x, y) until the next curve header.
	cc = contourc(gridx, gridy, img, [level, level]);
	
	%Begin at the first header
	cur_ind = 1;
	
	%The expected array is preallocated so we don't have to grow it every loop
	atoms_x = exp_n; atoms_y = exp_n;% atoms_r = exp_n;
	cur_atom = 0;
		
	while cur_ind <= length(cc(1,:))
		%Get the number of points in the curve and extract the points
		cur_accum = cc(2, cur_ind);
		cont_x = cc(1, cur_ind+1 : cur_ind+cur_accum)';
		cont_y = cc(2, cur_ind+1 : cur_ind+cur_accum)';
		
		%Least-squares circle fitting, adapted from circfit by Izhak Bucher, see license below
		p_ls = [cont_x, cont_y, ones(size(cont_x))] \ (-(cont_x.^2+cont_y.^2));
		p_x = -.5*p_ls(1); p_y = -.5*p_ls(2); %p_r = sqrt(-p_ls(3) + p_x^2 + p_y^2);
		
		%Update the return arrays
		cur_atom = cur_atom + 1;
		atoms_x(cur_atom) = p_x; atoms_y(cur_atom) = p_y; %atoms_r(cur_atom) = p_r;
		
		%Update the index for the next loop
		cur_ind = cur_ind + cur_accum + 1;
	end
	
	%Trim the arrays
	atoms_x = atoms_x(1:cur_atom); atoms_y = atoms_y(1:cur_atom); %atoms_r = atoms_r(1:cur_atom);
end


%circfit license
%Copyright (c) 1981, Izhak Bucher
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%* Redistributions of source code must retain the above copyright notice, this
%  list of conditions and the following disclaimer.
%
%* Redistributions in binary form must reproduce the above copyright notice,
%  this list of conditions and the following disclaimer in the documentation
%  and/or other materials provided with the distribution
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
%FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
