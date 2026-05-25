clc; clear; close all; %gm35_2_cd  >> OK both

%% 1. Select File
% [fileName, fileDir] = uigetfile('*.ffd', 'Select HFSS Pattern File');
% if isequal(fileName, 0), return; end
% filePath = fullfile(fileDir, fileName);

filePath = 'E:\Data\samples\antenna_pattern_sample_HFSS1.ffd';
filePath = 'E:\Data\samples\antenna_pattern_sample_HFSS.ffd';
% filePath = 'E:\Data\samples\antenna_pattern_sample_HFSS_2.ffd';

%% 2. Parse Grid Parameters (Handles floating-point angles safely)
fid = fopen(filePath, 'r');
if fid == -1, error('Failed to open file.'); end

th = num2cell(fscanf(fid, '%f', 3)); 
thetaStart = th{1}; thetaStop = th{2}; thetaCount = round(th{3}); 
thetaVec = linspace(thetaStart, thetaStop, thetaCount).';

ph = num2cell(fscanf(fid, '%f', 3));  
phiStart = ph{1}; phiStop = ph{2}; phiCount = round(ph{3}); 
phiVec = linspace(phiStart, phiStop, phiCount).';

nPts = thetaCount * phiCount;

% Read next non-empty line to find "frequencies" count safely
header_line = '';
while ischar(header_line) && isempty(strtrim(header_line))
    header_line = fgetl(fid);
end
fclose(fid);

% Identify file structure
if ischar(header_line) && contains(lower(header_line), 'frequencies')
    freqNum = sscanf(header_line, '%*s %d');
    isFreqDep = true;
    numHeaderLines = 3;
else
    freqNum = 0;
    isFreqDep = false;
    numHeaderLines = 2;
end

%% 3. Read Matrix Data
% Mfull = readmatrix(filePath, FileType="text", NumHeaderLines=numHeaderLines);
% Mfull = readmatrix(filePath, FileType="text", Delimiter = ["\t", " "], NumHeaderLines=0);
Mfull = readmatrix(filePath, FileType="text", Delimiter = ["\t", " "], NumHeaderLines=numHeaderLines);

% Extract frequencies (first column is NaN, second contains frequency)
freq_rows = Mfull(isnan(Mfull(:, 1)) & ~isnan(Mfull(:, 2)), :);
freqs = freq_rows(:, 2);

% Extract numeric field data
MData = Mfull(~isnan(Mfull(:, 1)), :);

if ~isFreqDep
    freqs = NaN; % Assign NaN for Frequency Independent files to unify logic
end
nBlk = numel(freqs);

assert(size(MData, 1) == nPts * nBlk, 'Parsed data points count (%d) does not match expected size (%d).', size(MData, 1), nPts * nBlk);

%% 4. Vectorized Split & Coordinate Table Construction
% Cell Method 2: Extremely simple, fast, and handles single/multi blocks
MCell = mat2cell(MData, repmat(nPts, nBlk, 1), size(MData, 2));

% Generate complete coordinate columns (Theta is outer loop, Phi is inner loop)
Theta_col = repelem(thetaVec, phiCount, 1);
Phi_col   = repmat(phiVec, thetaCount, 1);

% Convert each cell to a closed, sorted, 6-column table
MCellTables = cell(nBlk, 1);
for idx = 1:nBlk
    selected_field = MCell{idx};
    
    % Construct base table
    T = table(Theta_col, Phi_col, ...
        selected_field(:, 1), selected_field(:, 2), ...
        selected_field(:, 3), selected_field(:, 4), ...
        'VariableNames', {'Theta', 'Phi', 'E_theta_real', 'E_theta_imag', 'E_phi_real', 'E_phi_imag'});
    
    % Check if Phi is in the [-180, 180] range and requires shifting
    has_neg_phi = any(T.Phi < 0);
    if has_neg_phi
        if idx == 1
            fprintf('Detected Phi in [-180, 180] range. Shifting to [0, 360] range...\n');
        end
        % 1. Shift negative Phi values to [180, 360) range
        neg_idx = T.Phi < 0;
        T.Phi(neg_idx) = T.Phi(neg_idx) + 360;
        
        % 2. Deduplicate overlapping endpoints (e.g. -180 and 180 both mapping to 180)
        [~, unique_idx] = unique(T(:, {'Theta', 'Phi'}), 'rows', 'first');
        T = T(unique_idx, :);
        
        % 3. Sort primarily by Theta (ascending) and secondarily by Phi (ascending)
        T = sortrows(T, {'Theta', 'Phi'});
    end
    
    % 4. Close the sphere/seam by copying Phi = 0 to Phi = 360
    if ~any(T.Phi == 360)
        if idx == 1 && has_neg_phi
            fprintf('Closing the sphere seam: copying Phi = 0 to Phi = 360...\n');
        end
        rows_phi0 = T(T.Phi == 0, :);
        rows_phi360 = rows_phi0;
        rows_phi360.Phi(:) = 360;
        T = [T; rows_phi360];
        T = sortrows(T, {'Theta', 'Phi'});
    end
    
    MCellTables{idx} = T;
end

%% 5. Struct Method 1: Bind frequencies and closed-sphere tables together
MStruct = struct('blockID', num2cell((1:nBlk).'), 'freq', num2cell(freqs), 'data', MCellTables);

%% 6. Access Struct Data and Plot (e.g., First Frequency Block)
plot_idx = 1; 
selected_table = MStruct(plot_idx).data;
selected_freq  = MStruct(plot_idx).freq;

% Dynamically extract grid dimensions directly from the processed table
thetaVec_new = unique(selected_table.Theta);
phiVec_new   = unique(selected_table.Phi);
thetaCount_new = numel(thetaVec_new);
phiCount_new   = numel(phiVec_new);

% Reconstruct complex E-fields directly from the struct table
E_theta_vec = selected_table.E_theta_real + 1i * selected_table.E_theta_imag;
E_phi_vec   = selected_table.E_phi_real   + 1i * selected_table.E_phi_imag;

% Reshape 1D data to 2D grids using the dynamic dimensions
E_theta_grid = reshape(E_theta_vec, [phiCount_new, thetaCount_new]).';
E_phi_grid   = reshape(E_phi_vec, [phiCount_new, thetaCount_new]).';

[Phi_grid, Theta_grid] = meshgrid(phiVec_new, thetaVec_new);

% Calculate Total Gain
E_total = sqrt(abs(E_theta_grid).^2 + abs(E_phi_grid).^2);
Gain_dB = 20 * log10(E_total);
max_gain = max(Gain_dB(:));
Gain_clamped = max(Gain_dB, max_gain - 40); % Cap dynamic range to 40 dB

% Calculate Circular Polarization & Axial Ratio
E_rhcp = (E_theta_grid + 1i * E_phi_grid) / sqrt(2);
E_lhcp = (E_theta_grid - 1i * E_phi_grid) / sqrt(2);
AR_linear = (abs(E_rhcp) + abs(E_lhcp)) ./ max(abs(abs(E_rhcp) - abs(E_lhcp)), 1e-12);
AR_dB = min(20 * log10(AR_linear), 40); % Cap display at 40 dB

%% 7. Render Layout Dashboard
freq_label = iff(isFreqDep, sprintf(' @ %.2f GHz', selected_freq/1e9), ' (Freq. Independent)');
figure('Name', 'Antenna Performance Dashboard', 'Color', 'w', 'Position', [100, 100, 1100, 750]);

% 2D Contour: Total Gain
subplot(2, 2, 1);
contourf(Phi_grid, Theta_grid, Gain_clamped, 25, 'LineColor', 'none');
colorbar; colormap(gca, 'jet'); grid on;
xlabel('\Phi (deg)'); ylabel('\Theta (deg)');
title(['Total Gain (dB)', freq_label]);
set(gca,'YDir','reverse');

% 2D Contour: Axial Ratio
subplot(2, 2, 2);
contourf(Phi_grid, Theta_grid, AR_dB, 0:0.5:10, 'LineColor', 'none');
colorbar; colormap(gca, 'jet'); clim([0, 10]); grid on;
xlabel('\Phi (deg)'); ylabel('\Theta (deg)');
title(['Axial Ratio (dB)', freq_label]);
set(gca,'YDir','reverse');

% 3D Coordinate Mapping
R_3D = Gain_clamped - min(Gain_clamped(:)); 
X = R_3D .* sind(Theta_grid) .* cosd(Phi_grid);
Y = R_3D .* sind(Theta_grid) .* sind(Phi_grid);
Z = R_3D .* cosd(Theta_grid);

% 3D Spherical Plot: Total Gain
subplot(2, 2, 3);
surf(X, Y, Z, Gain_clamped, 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
colorbar; colormap(gca, 'jet'); axis equal; grid on; view(135, 30); camlight;
xlabel('X (dB)'); ylabel('Y (dB)'); zlabel('Z (dB)');
title('3D Total Gain (dB)');

% 3D Spherical Plot: Axial Ratio mapped to Gain Structure
subplot(2, 2, 4);
surf(X, Y, Z, AR_dB, 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
colorbar; colormap(gca, 'jet'); clim([0, 10]); axis equal; grid on; view(135, 30); camlight;
xlabel('X (dB)'); ylabel('Y (dB)'); zlabel('Z (dB)');
title('3D Pattern Colored by Axial Ratio (dB)');

% Simple inline conditional helper
function out = iff(condition, true_val, false_val)
    if condition, out = true_val; else, out = false_val; end
end