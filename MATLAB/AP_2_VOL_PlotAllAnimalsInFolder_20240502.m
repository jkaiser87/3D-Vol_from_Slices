% This script summarized all animals within a given folder and plots them
% into one 3D brain (flipped if necessary)

% open folder that contains either "*_variables.mat" files or subfolders
% with animals that have been processed 

ExperimentName = 'Exp'; %as prefix for files
fliptoright = 0; %flips the volume to the right hemisphere if wanted, otherwise put to 0 to keep it at the side it was traced at
fliptoleft = 1; %flips the volume to the right hemisphere if wanted, otherwise put to 0 to keep it at the side it was traced at

%%
% Define base directory and find files
baseDir = pwd;
files = dir(fullfile(baseDir, '**', '*_variables.mat'));

% Prepare color scheme
baseColors = {[1, 0, 0], [0, 1, 0], [0, 0, 1]};  % Red, Green, Blue for channels C1, C2, etc.

%colorAdjustmentFactor = linspace(0.5, 1, numel(files)); % Generate lightness adjustments
hueIncrement = linspace(-0.05, 0.05, numel(files));  % Small increments around the base hue
saturationFactor = linspace(0.6, 1, numel(files));  % Ensure this stays within 0 to 1 after adjustment
alpha = 0.1; %how transparent each volume should be

% Initialize figures
figure1 = figure('Name', 'All Animals Combined', 'Color', 'w');
fig1ax = axes;
set(fig1ax,'ZDir','reverse');

figure2 = figure('Name', 'Individual Animals', 'Color', 'w');
fig2ax = axes;
set(fig2ax,'ZDir','reverse');
numAnimals = length(files);
numRows = ceil(sqrt(numAnimals));
numCols = ceil(numAnimals / numRows);

midline=575;

% Prepare the brain outline for both figures
allen_atlas_path = fullfile(userpath, 'AP_histology', 'allenAtlas');
av = readNPY(fullfile(allen_atlas_path, 'annotation_volume_10um_by_index.npy'));
% Define mesh parameters
slice_spacing = 5;
reduced_av = av(1:slice_spacing:end, 1:slice_spacing:end, 1:slice_spacing:end);
brain_volume = bwmorph3(bwmorph3(reduced_av > 1,'majority'), 'majority');
brain_outline_patchdata = isosurface(permute(brain_volume,[3,1,2]), 0.5);

% Plot brain in figure 1 (combined plot)
figure(figure1);
patch('Vertices', brain_outline_patchdata.vertices * slice_spacing, ...
    'Faces', brain_outline_patchdata.faces, ...
    'FaceColor', [0.7,0.7,0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
hold on;
view([90,90]);
axis('vis3d','equal','off','manual');
axis tight;

legendHandles = [];
legends = [];
% Number of rows and columns for the subplot
numRows = 1;
numCols = numAnimals;

% Create a figure for the individual subplots
figure2 = figure;

% Set the subplot layout with minimal space between panels
tiledlayout(numRows, numCols, 'TileSpacing', 'compact', 'Padding', 'compact');

% Process each animal
for idx = 1:numAnimals
    fullPath = fullfile(files(idx).folder, files(idx).name);
    [~, animalName, ~] = fileparts(files(idx).name);
    animalName = erase(animalName, '_variables');

    % Load data
    load(fullPath, 'inj_vol');

    % Plot in individual subplot
    nexttile;
    patch('Vertices', brain_outline_patchdata.vertices * slice_spacing, ...
    'Faces', brain_outline_patchdata.faces, ...
    'FaceColor', [0.7,0.7,0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    hold on;
    view([90,90]);
    axis('vis3d','equal','off','manual');
    axis tight;

    % Plot each channel for both figures
    for ch = 1:length(inj_vol)
        channelData = inj_vol(ch);
        channelLabel = channelData.channel{1};  % Assuming channel names are stored as strings in cell arrays
        vertices = channelData.smoothedVertices;

        if (fliptoright == 1)
            y_coords = vertices(:,2);
            flip_indices = y_coords < midline;
            vertices(flip_indices, 2) = 2 * midline - y_coords(flip_indices);
        end

        if (fliptoleft == 1)
            y_coords = vertices(:,2);
            flip_indices = y_coords > midline;
            vertices(flip_indices, 2) = 2 * midline - y_coords(flip_indices);
        end

        % Determine the base color for the channel
        colorIndex = mod(ch-1, numel(baseColors)) + 1;
        baseColor = baseColors{colorIndex};
        
        % Adjust color
        hsvColor = rgb2hsv(baseColor);
        hsvColor(1) = mod(hsvColor(1) + hueIncrement(idx), 1);  % Adjust hue, wrap around using mod
        hsvColor(2) = min(1, max(0, hsvColor(2) * saturationFactor(idx)));  % Adjust saturation, ensuring it stays within [0,1]
        %hsvColor(3) = hsvColor(3) * colorAdjustmentFactor(idx);  % Adjust brightness
        adjustedColor = hsv2rgb(hsvColor);

        % Plot in the combined figure
        figure(figure1);
        h = patch('Vertices', vertices, 'Faces', channelData.k1, 'FaceColor', adjustedColor, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        % Plot in the individual subplot
        figure(figure2);
        patch('Vertices', vertices, 'Faces', channelData.k1, 'FaceColor', adjustedColor, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
        % Collect legend names
        legendHandles(end+1) = h;  % Add the handle to the array
        legends{end+1} = sprintf('%s - %s', channelLabel, animalName);
    end

    % Add title to the subplot
    title(animalName);
end

% Finalize combined figure settings and saving
figure(figure1);
xlabel('AP');
ylabel('ML');
zlabel('DV');

% Sorting the legend
legendInfo = struct('handle', num2cell(legendHandles), 'label', legends);
[~, idx] = sort({legendInfo.label});
sortedLegendInfo = legendInfo(idx);
sortedHandles = [sortedLegendInfo.handle];
sortedLabels = {sortedLegendInfo.label};

% Adding the sorted legend
%legend(sortedHandles, sortedLabels, 'Location', 'northeastoutside', 'NumColumns', 2);
view(90, 90);
axis tight;
rotate3d on;
%% uncomment if you want to also capture a video turn (need to install CaptureFigVid)

%OptionZ.FrameRate=60;OptionZ.Duration=5.5;OptionZ.Periodic=true;
%CaptureFigVid([90,90;-50,45;90,90], 'OverviewAllAnimals',OptionZ)

%%

% Save the current figure
savename=fullfile(baseDir,[ExperimentName,'_3DPlot_All_Animals_Combined']);
savefig(gcf,savename);
print([savename,'.png'], '-dpng');
print('-dpng', '-r300', fullfile(baseDir,[ExperimentName,'_3DPlot_All_Animals_Combined_HighRes.png']));  % Save as High-Resolution PNG

disp(['Figure 1 saved as ',savename,'.png and .m']);

view([-50,45]);
savename=fullfile(baseDir,[ExperimentName,'_3DPlot_All_Animals_Combined_angled']);
savefig(gcf,savename);
print([savename,'.png'], '-dpng');
print('-dpng', '-r300', fullfile(baseDir,[ExperimentName,'_3DPlot_All_Animals_Combined_angled_HighRes.png']));  % Save as High-Resolution PNG


% Finalize subplot figure settings
figure(figure2);
view(90, 90);
axis tight;
rotate3d on;

savename=fullfile(baseDir,[ExperimentName,'_3DPlot_Individual_Animals']);
savefig(gcf,savename);
print([savename,'.png'], '-dpng');
disp(['Figure 2 saved as ',savename,'.png and .m']);
print('-dpng', '-r300', fullfile(baseDir, 'All_Animals_Combined_HighRes.png'));  % Save as High-Resolution PNG

