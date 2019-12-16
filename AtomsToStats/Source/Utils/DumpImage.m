%This saves the exact index of pixel array img to file fname using the colormap cmap. NaN entries are colored black
%Gabriel Martine
%Last updated: Jul 27 2017

function DumpImage(img, fname, cmap, vmin, vmax, res)
	%We must first convert the MxN matrix img to an MxNx3 array of RGB data
	%This is done by clipping the array to res colors including vmin and vmax. The edges
	%are set to be equidistant to levels of the color ladder so that the values vmin+eps
	%and vmax-eps are correctly mapped to vmin and vmax instead of vmin+step and vmax, as
	%in Matlab's default coloring routines
	cmap = eval(strcat(cmap, '(', num2str(int16(res)), ')'));
	cmap(end+1, :) = [0.0, 0.0, 0.0];
	edges = vmin + (0.5 + (0:res-2))*(vmax-vmin)/(res-1);
	
	inv = isnan(img);
	ind = img;
	for nn = 1:res-1
		msk = (img <= edges(nn)); img(msk) = NaN;
		ind(msk) = nn;
	end
	ind(~isnan(img)) = res;
	ind(inv) = res+1;
	
	imwrite(ind2rgb(flipud(ind), cmap), fname)
end
