%Creates a time index array to be used in the evolution for energy calculation and phase saving steps
%Gabriel Martine
%Last updated: Dec 17 2016

function time_ind = SaveTimes(begin_ind, end_ind, snap)
	if snap <= 0
		%A non-positive value for snap ignores the step
		time_ind = 0;
	else
		if begin_ind == 0
			%Choose linear spacing in linear time if the starting index is 0
			time_ind = round(linspace(0, end_ind, snap+1));
			time_ind(1) = [];
		else
			%Choose linear spacing in logarithmic time otherwise
			time_ind = round(begin_ind*(end_ind/begin_ind).^((1:snap)/snap));
		end
		
		%All indices must be larger than begin_ind since saving checks happen after begin_ind
		time_ind(time_ind <= begin_ind) = begin_ind + 1;
		
		%Remove redundant time indices
		time_ind = unique(time_ind);
	end
end
