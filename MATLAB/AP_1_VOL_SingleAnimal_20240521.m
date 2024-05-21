%% Volume analysis

% Julia Kaiser, April 2024
% this script is a little helper to process injection volumes traced in
% FIJI and plot them into 3D space using AP_histology written by petersaj
% https://github.com/petersaj/AP_histology
% install AP-histology according to official guidelines

% Preparations before running this script:

% --- File naming convention: EXP-ANIMAL_...._s001.tif etc
% ----- use Bulk rename utility to fix if necessary (- and _ are
% important!) - DO THIS BEFORE FIJI!!!

% --- run through FIJI toolset to create slices, rotate & flip, and trace injection volumes
% ----- make sure that there is only 1 animal in folder 
% should contain subfolder VOL/ and *.tif slices in main folder (open)



% adapt these settings:
rerun_histology = 1; %redo AP_Histology? can be 0 if already done
channelsToProcess = {'C1', 'C2'}; % Define which channels to process (should have csv files in subfolder)
colors={'red','green'}; %colors to plot the channels by
structure_names = {'Somatomotor areas', 'Somatosensory areas', ...
    'Visual areas', 'Auditory areas'}; %structures to plot into the brain (light grey)
%structure_names = {}; %uncomment this line if you dont want to plot any brain structures


%%
%----- NO CHANGES HERE
parentFolder = pwd;
cd(parentFolder);

subfolderName = 'OUT';

% Full path to the new subfolder
outputFolderPath = fullfile(parentFolder, subfolderName);

% Check if the subfolder already exists
if ~exist(outputFolderPath, 'dir')
    % The subfolder does not exist, so create it
    mkdir(outputFolderPath);
end

figPath = fullfile(outputFolderPath, 'FIG');
if ~exist(figPath, 'dir')
    % The subfolder does not exist, so create it
    mkdir(figPath);
end


csvsavePath = fullfile(outputFolderPath, 'CSV');
if ~exist(csvsavePath, 'dir')
    % The subfolder does not exist, so create it
    mkdir(csvsavePath);
end

% Calculate %ap for min and max
z_min = 8; % Min limit from X-axis
z_max = 1320; % Max limit from X-axis
y_min = 58;
y_max = 1093;

midline = 575;

%%
% run through the following pipeline
% DO NOT DOWNSAMPLE
% DO NOT PROCESS IMAGES (rotate/flip/...)

if (rerun_histology==1)
    AP_histology
    disp('Go through all steps in AP_histology window, then continue by pressing any key inside the terminal. \n To stop, press Ctrl + C');
    pause;
end

%% data processing to convert volumes and plot

% --- GATHER COORDINATES INTO STRUCT
% Get a list of all CSV files in the folder
folderPath = fullfile(parentFolder,'VOL/CSV'); %filepath to csv files

[parentFolderPath, parentFolderName, ~] = fileparts(parentFolder);
%[~, parentFolderName, ~] = fileparts(parentFolderPath);
disp(['-------- Processing ',parentFolderName]); %second folder up to Slices

% Get a list of all TIFF files in the folder
tifFiles = dir(fullfile(parentFolder, '*.tif'));

% Initialize a structure to store points, dynamically creating fields for each channel
points = struct;

% Loop through each TIFF file
for i = 1:length(tifFiles)
    baseFileName = tifFiles(i).name;
    [~, name, ~] = fileparts(baseFileName);

    % Process each specified channel
    for j = 1:length(channelsToProcess)
        channel = channelsToProcess{j};
        if ~isfield(points, channel)
            points.(channel) = struct('name', {}, 'X', {}, 'Y', {});
        end
        csvPath = fullfile(folderPath, strcat(channel, '_', name, '.csv'));

        points.(channel)(i).name = name;
        points.(channel)(i).X = [];
        points.(channel)(i).Y = [];
        % Check if the CSV file exists
        if exist(csvPath, 'file')
            dataTable = readtable(csvPath);
            points.(channel)(i).X = dataTable.X;
            points.(channel)(i).Y = dataTable.Y;
        end
    end
end

disp('Files read into dataTable');

savemat=fullfile(outputFolderPath,[parentFolderName,'_variables.mat']);
save(savemat,'points'); 

%%

% --- CONVERT POINTS TO CCF
% Load histology/CCF alignment
ccf_slice_fn = fullfile(outputFolderPath,'histology_ccf.mat');
load(ccf_slice_fn);
ccf_alignment_fn = fullfile(outputFolderPath,'atlas2histology_tform.mat');
load(ccf_alignment_fn);

ccf_points = struct;

for j = 1:length(channelsToProcess)
    channel = channelsToProcess{j};
    dataSets = points.(channel); %subset to this specific channel
    ccf_points.(channel) = cell(length(dataSets), 1);

    for i = 1:length(dataSets)
        histology_points = [dataSets(i).X, dataSets(i).Y];

        if ~isempty(histology_points)
            % Transform histology to atlas slice
            tform = affine2d;
            tform.T = atlas2histology_tform{i}; % Adjust index i according to your data
            tform = invert(tform);

            % Transform and round to nearest index
            [histology_points_atlas_x,histology_points_atlas_y] = ...
                transformPointsForward(tform, ...
                histology_points(:,1), ...
                histology_points(:,2));

            histology_points_atlas_x = round(histology_points_atlas_x);
            histology_points_atlas_y = round(histology_points_atlas_y);

            [M, N] = size(histology_ccf(i).av_slices);
            outOfRangeY = histology_points_atlas_y < 1 | histology_points_atlas_y > M;
            outOfRangeX = histology_points_atlas_x < 1 | histology_points_atlas_x > N;

            if any(outOfRangeY)
                fprintf('Y indices out of range: Min = %d, Max = %d, Valid Range = [1, %d]\n', min(histology_points_atlas_y(outOfRangeY)), max(histology_points_atlas_y(outOfRangeY)), M);
            end

            if any(outOfRangeX)
                fprintf('X indices out of range: Min = %d, Max = %d, Valid Range = [1, %d]\n', min(histology_points_atlas_x(outOfRangeX)), max(histology_points_atlas_x(outOfRangeX)), N);
            end


            probe_points_atlas_idx = sub2ind(size(histology_ccf(i).av_slices), histology_points_atlas_y,histology_points_atlas_x);

            % Get CCF coordinates for histology coordinates (CCF in AP,DV,ML)
            ccf_points.(channel){i} = ...
                [histology_ccf(i).plane_ap(probe_points_atlas_idx), ...
                histology_ccf(i).plane_dv(probe_points_atlas_idx), ...
                histology_ccf(i).plane_ml(probe_points_atlas_idx)];

        end
    end
end

disp('Coordinates transferred into CCF space');
save(savemat,'ccf_points','-append'); 

% Outputs
% ccf_points = cell array with CCF coordinates corresponding to histology_points (note: in native CCF order [AP/DV/ML])

%% plot

allen_atlas_path = fullfile(userpath,'AP_histology\allenAtlas');
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

figure('Color','w');
ccf_3d_axes = axes;
set(ccf_3d_axes,'ZDir','reverse');
hold(ccf_3d_axes,'on');
axis(ccf_3d_axes,'vis3d','equal','off','manual');
view([90,90]);
axis tight;
rotate3d on;

% plot brain
% Define mesh parameters
slice_spacing = 5; % Adjust based on your resolution needs
reduced_av = av(1:slice_spacing:end, 1:slice_spacing:end, 1:slice_spacing:end);

% Generate a binary volume where non-zero voxels are considered part of the brain
brain_volume = bwmorph3(bwmorph3(reduced_av > 1,'majority'), 'majority');

% Generate the mesh using isosurface
brain_outline_patchdata = isosurface(permute(brain_volume,[3,1,2]), 0.5);

% Plot the mesh as a patch on the 3D axes
patch(ccf_3d_axes, ...
    'Vertices', brain_outline_patchdata.vertices * slice_spacing, ...
    'Faces', brain_outline_patchdata.faces, ...
    'FaceColor', [0.7,0.7,0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
hold on; % Hold on to plot multiple categories

% Define distinct grey tones
num_structures = length(structure_names);
grey_values = linspace(0.2, 0.8, num_structures);  % Create a gradient from light grey to dark grey

% Initialize cell arrays for legend handling
legends = cell(num_structures, 1);
patch_handles = zeros(num_structures, 1);

% Plot the structures with distinct grey tones and black edges

struct_vol = struct;

for i = 1:num_structures
    structure_name = structure_names{i};
    structure_index = find(strcmpi(st.safe_name, structure_name), 1);

    if isempty(structure_index)
        warning(['Structure not found: ', structure_name]);
        continue;
    end

    target_path = st.structure_id_path{structure_index}; % Get the specific path component outside of cellfun
    plot_ccf_idx = find(cellfun(@(x) contains(x, target_path), st.structure_id_path)); % Use the predefined path
    plot_ccf_volume = ismember(reduced_av, plot_ccf_idx);
    structure_3d = isosurface(permute(plot_ccf_volume, [3, 1, 2]), 0);

    if ~isempty(structure_3d.vertices)
        struct_vol(i).name = {structure_name};
        struct_vol(i).struct = structure_3d;
        h = patch(ccf_3d_axes, 'Vertices', structure_3d.vertices * slice_spacing, 'Faces', structure_3d.faces, ...
            'FaceColor', repmat(grey_values(i), 1, 3), 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        patch_handles(i) = h;
        legends{i} = st.safe_name{structure_index};
    else
        warning(['No visible structure for plotting: ', st.safe_name{structure_index}]);
    end
end


inj_vol = struct;
% PLOT VOLUMES
for j = 1:length(channelsToProcess)
    channel = channelsToProcess{j};
    dataSets = ccf_points.(channel); %subset to this specific channel
    ccf_points_cat = round(cell2mat(dataSets)); %z, x, y
    ccf_points_cat_ord = [ccf_points_cat(:,1),ccf_points_cat(:,3),ccf_points_cat(:,2)];
    ccf_points_cat_ord = interpolatePoints(ccf_points_cat_ord, 50);  % Increase '50' as needed for more density
    [k1, vol] = boundary(ccf_points_cat_ord);
    smoothedVertices = laplacianSmooth(ccf_points_cat_ord, k1, 0.1, 5); % Apply the smoothing
    inj_vol(j).channel = {channel};
    inj_vol(j).smoothedVertices = smoothedVertices;
    inj_vol(j).k1 = k1;
    inj_vol(j).ccf_points_cat_ord = ccf_points_cat_ord;
    h = patch('Vertices', smoothedVertices, 'Faces', k1, 'FaceColor', colors{j}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch_handles(end+1) = h;
    legends{end+1} = colors{j};
end


% Create the legend
h = legend(patch_handles(patch_handles > 0), legends{patch_handles > 0}, 'Location', 'BestOutside');
hold on
view(90, 90);

% Set axis labels and adjust plot properties as needed
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

% Save the current figure
savename=fullfile(subfolderName,'FIG',[parentFolderName,'_3DPlot_Volume']);
savefig(gcf,savename);
% Save the current figure to a PNG file
print([savename,'.png'], '-dpng');
disp(['Figure saved as ',savename,'.png and .m']);

save(savemat,'inj_vol','struct_vol','-append'); 


% For high-resolution, you can specify the resolution in DPI (dots per inch)
%print([savename,'_highres.png'], '-dpng', '-r300'); % 300 DPI

%% calculate mid and SD

FinalData = array2table(zeros(0,8), 'VariableNames', {'ap', 'dv', 'ml', 'coord', 'Channel', 'apperc', 'mlperc', 'hemisphere'});
currentData = array2table(zeros(0,8), 'VariableNames', {'ap', 'dv', 'ml', 'coord', 'Channel', 'apperc', 'mlperc', 'hemisphere'});

% Loop through channels
for j = 1:length(channelsToProcess)
    channel = channelsToProcess{j};
    dataSets = ccf_points.(channel); %subset to this specific channel
    ccf_points_cat = round(cell2mat(dataSets));

    % Initialize min and max coordinates
    min_coords = [inf, inf, inf];
    max_coords = [-inf, -inf, -inf];

    % Loop through all vertices to find min and max coordinates
    for i = 1:size(ccf_points_cat, 1)
        % Update min and max coordinates along each axis
        min_coords = min(min_coords, ccf_points_cat(i, :));
        max_coords = max(max_coords, ccf_points_cat(i, :));
    end

    % Create a table for min and max coordinates
    minmax = vertcat(min_coords, max_coords);
    currentData = array2table(minmax, 'VariableNames', {'ap', 'dv', 'ml'});

    % Add coord, Channel, apperc, mlperc, and hemisphere columns
    currentData.coord = {'min'; 'max'};
    currentData.Channel = repmat({channel}, 2, 1);
    currentData.apperc = NaN(2, 1);
    currentData.mlperc = NaN(2, 1);
    currentData.hemisphere = repmat({'L'}, 2, 1);

    % Calculate % in z values
    normalized_z = ((currentData.ap - z_min) / (z_max - z_min)) * 100;
    currentData.apperc = normalized_z;

    % Calculate % to midline
    normalized_y = ((currentData.ml - y_min) / (y_max - y_min)) * 100;
    currentData.mlperc = normalized_y;

    % Determine hemisphere
    currentData.hemisphere(currentData.ml > midline) = {'R'};

    % Append currentData to data
    FinalData = [FinalData; currentData];
end


disp('Min/Max of all axes calculated');
disp('% in ap and ml calculated');

savename=fullfile(csvsavePath,[parentFolderName,'_MinMaxCoordsPerc']);
writetable(FinalData,[savename,'.csv']);
disp(['Table saved as ',savename,'.csv']);

%%
% requirements to run this:
% --- install the Curve Fitting Toolbox (through Edwin)
% --- the natsortfile add-on
% --- install https://github.com/petersaj/AP_histology and clone it into a
% MATLAB path folder (eg Documents > MATLAB, or somewhere else and add it
% to your folders) - this includes all dependencies, eg npy-matlab

disp('Done.');