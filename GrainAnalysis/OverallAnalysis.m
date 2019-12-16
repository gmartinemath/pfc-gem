%%%%%%%%%%%%%% This code reads extracted grain distributions and outputs graphs for the
%%%%%%%%%%%%%% statistical distributions
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: Apr 10 2019

clear all;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options
area_to_diameter = true;
pfc_min_coord = 2; pfc_max_coord = 14;
out_folder = '~/Downloads/PFC_Test/Stats_Run/';
legend_labels = {'PFC 9', 'PFC 10',};
pfc_color = {'b', 'r'};
pfc_mark = 's';
thickness = 1.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input files
PFC_IN_FOLDER = '~/Downloads/PFC_Test/ATS_Run/';
PFC_SUBDIRS = 1:1;
PFC_INDICES = [9, 10];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create the figures
figure(1); hold on; box on;
if area_to_diameter
	title(''); ylabel('Probability Density'); xlabel('Reduced Equivalent Diameter'); xlim([0, 3]);
else
	title(''); ylabel('Probability Density'); xlabel('Reduced Area'); xlim([0, 5]);
end

figure(2); hold on; box on; xlim([0, 3]); title(''); ylabel('Probability Density'); xlabel('Reduced Perimeter');
figure(3); hold on; box on; xlim([0.5, 1]); title(''); ylabel('Probability Density'); xlabel('Isoperimetric Ratio');
figure(4); hold on; box on; xlim([0, 0.2]); title(''); ylabel('Probability Density'); xlabel('Convex Hull Ratio');

figure(5); hold on; box on; xlim([1, pfc_max_coord+1]); title(''); ylabel('Number Fraction'); xlabel('Number of Sides');

figure(6); hold on; box on; xlim([1, pfc_max_coord+1]); title(''); ylabel('Average Area of Side Class');
xlabel('Number of Sides');

figure(7); hold on; box on; xlim([1, pfc_max_coord+1]); title(''); ylabel('Area fraction of Side Class');
xlabel('Number of Sides');

figure(8); hold on; box on; xlim([1, pfc_max_coord+1]); title('');
ylabel('Average Side Class of Neighbors'); xlabel('Number of Sides');

%figure(9); hold on; box on; title('Lognormal Mu Fit Parameter Comparison'); xlabel('PFC Grain Number'); ylabel('Mu');
%figure(10); hold on; box on; title('Lognormal Sigma Fit Parameter Comparison'); xlabel('PFC Grain Number'); ylabel('Sigma');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PFC data
logn_x = []; logn_mu = []; logn_err_mu = []; logn_sigma = []; logn_err_sigma = [];
for pfc_ind = 1:length(PFC_INDICES)
	%Accumulate all statistics
	stats = 0;
	for pfc_subdirs = 1:length(PFC_SUBDIRS)
	stats = MergeGrainStats(stats, ...
		GetGrainStats(sprintf('%sG%d/Data_%d.mat', PFC_IN_FOLDER, PFC_SUBDIRS(pfc_subdirs), PFC_INDICES(pfc_ind))));
	end
	PlotGrainStats(stats, area_to_diameter, pfc_min_coord, pfc_max_coord, pfc_color{pfc_ind}, thickness, pfc_mark);
	
	%Get PFC's lognormal parameters
	%[mu, sigma, err_mu, err_sigma] = LognormalParams(sqrt(stats.area)/mean(sqrt(stats.area)));
	%logn_x(end+1) = length(stats.area)/length(PFC_SUBDIRS);
	%logn_mu(end+1) = mu; logn_err_mu(end+1) = err_mu;
	%logn_sigma(end+1) = sigma; logn_err_sigma(end+1) = err_sigma;
	%fprintf('MU: %.2f (%.2f)     SIGMA: %.2f (%.2f)', mu, err_mu, sigma, err_sigma);
	%exp(mu+0.5*sigma^2)
end	

%Plot the lognormal fit parameters (reduced area only)
%figure(9); errorbar(logn_x, logn_mu, logn_err_mu, 'color', 'r', 'linewidth', thickness, 'LineStyle','none');
%figure(10); errorbar(logn_x, logn_sigma, logn_err_sigma, 'color', 'r', 'linewidth', thickness, 'LineStyle','none');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add the legend, use an appropriate scale and save
mkdir(out_folder);
save_str = {'ReducedArea', 'ReducedPerim', 'IPR', 'CHR', 'NumberOfSides', 'AAvsNS', 'AFSCvsNS', 'ANSvsNS', 'LogMu', 'LogSigma'};
if area_to_diameter
	save_str{1} = 'ReducedEquivalentDiameter';
end
for fig_id = 1:8%10
	figure(fig_id);
	legend_location = 'northeast';
	if ismember(fig_id, [3, 6])
		legend_location = 'northwest';
	elseif fig_id == 9
		legend_location = 'southeast';
	end
	legend(legend_labels, 'location', legend_location);
	
	%Increase the fontsize, remove the whitespace and save
	set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16)
 	ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset;
 	ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2), outerpos(3)-ti(1)-ti(3), outerpos(4)-ti(2)-ti(4)];
	saveas(fig_id, sprintf('%s%s.png', out_folder, save_str{fig_id}));
end
