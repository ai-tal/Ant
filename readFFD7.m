% Select an HFSS .ffd file %gm35_1 >> OK both
[fileName, fileDir] = uigetfile('*.ffd', 'Select an HFSS .ffd file');
if isequal(fileName,0) || isequal(fileDir,0), return; end
filePath = fullfile(fileDir, fileName);

% 1. Header Parsing (Robust fscanf to handle decimal angles)
fid = fopen(filePath, 'r');
th = num2cell(fscanf(fid, '%f', 3)); % start, stop, count
thetaStart = th{1}; thetaStop = th{2}; thetaCount = th{3};
thetaVec = linspace(th{:});

ph = num2cell(fscanf(fid, '%f', 3)); % start, stop, count
phiStart = ph{1}; phiStop = ph{2}; phiCount = ph{3};
phiVec = linspace(ph{:});

% Safely extract the optional Frequencies declaration line
C1 = textscan(fid, '%s %d', 1);
fclose(fid);

freqNum = 0;
if ~isempty(C1) && length(C1) >= 2 && ~isempty(C1{2})
    freqNum = C1{2};
end

isFreqDep = freqNum > 0;
numHeaderLines = 3;
if ~isFreqDep
    numHeaderLines = 2;
end

% 2. Bulk File Reading
Mfull = readmatrix(filePath, FileType="text", Delimiter = ["\t", " "], NumHeaderLines=numHeaderLines);

% 3. Extract Frequencies and Filter Fields
if isFreqDep
    freqs = Mfull(isnan(Mfull(:,1)), 2);
    nBlk = length(freqs);
else
    freqs = NaN;
    nBlk = 1;
end

% Extract only numeric field rows (skips rows starting with NaN "frequency")
MData = Mfull(~isnan(Mfull(:,1)), :);

% 4. Multi-Dimensional Vectorized Setup
MData_split = reshape(MData, [phiCount, thetaCount, nBlk, 4]);

% =========================================================================
% 5. Check & Rearrange Phi Coordinates to 0-360 (With Seam Closure)
% =========================================================================
if phiStart < 0
    % A. Map to [0, 360) and use unique() to resolve duplicate 180s
    phi_mod = mod(phiVec(:), 360);
    [phiVec, sort_idx] = unique(phi_mod); % unique sorts and filters duplicates
    
    % B. Apply sorting and duplicate reduction to the 4D matrix
    MData_split = MData_split(sort_idx, :, :, :);
    
    % C. Close the Seam: copy/append Phi=0 data to Phi=360 at the end
    phiVec = [phiVec; 360];
    phi_0_slice = MData_split(1, :, :, :);
    MData_split = cat(1, MData_split, phi_0_slice);
    
    % D. Update boundaries and counts
    phiCount = length(phiVec);
    phiStart = phiVec(1);
    phiStop  = phiVec(end);
else
    % If already positive, keep existing order
    sort_idx = 1:phiCount;
end

% 6. Flatten MData_split back to 2D in a single vectorized call
MData_clean = reshape(MData_split, [phiCount * thetaCount * nBlk, 4]);

% 7. Mapped Coordinate Columns
nPts = phiCount * thetaCount;
Theta_col = repelem(thetaVec(:), phiCount);
Phi_col   = repmat(phiVec(:), thetaCount, 1);

% =========================================================================
% 8. Vectorized Block Splitting & Struct Packing (Your Optimized Methods)
% =========================================================================
% Split the 2D matrix into cells using Cell Method 2
MCell = mat2cell(MData_clean, repmat(nPts, nBlk, 1), 4);

% Pre-generate the 6-column tables for each cell block
Tables = cell(nBlk, 1);
for f = 1:nBlk
    Tables{f} = table(Theta_col, Phi_col, ...
        MCell{f}(:, 1), MCell{f}(:, 2), ...
        MCell{f}(:, 3), MCell{f}(:, 4), ...
        'VariableNames', {'Theta', 'Phi', 'E_theta_real', 'E_theta_imag', 'E_phi_real', 'E_phi_imag'});
end

% Assemble the final structured array using Struct Method 1
MStruct = struct('blockID', num2cell((1:nBlk).'), ...
                 'freq',    num2cell(freqs), ...
                 'data',    MCell, ...
                 'table',   Tables);

% Print a preview of the structured array
fprintf('\nConstructed Struct Array (nBlk: %d, points per block: %d):\n', nBlk, nPts);
disp(MStruct(1));

% =========================================================================
% 9. Visualization (Reconstruct 2D arrays directly from the Table)
% =========================================================================
f_idx = 1; % Index of the table inside MStruct to visualize
T = MStruct(f_idx).table;

% Map columns back to standard 2D matrices [thetaCount, phiCount]
E_theta_real = reshape(T.E_theta_real, [phiCount, thetaCount]).';
E_theta_imag = reshape(T.E_theta_imag, [phiCount, thetaCount]).';
E_phi_real   = reshape(T.E_phi_real,   [phiCount, thetaCount]).';
E_phi_imag   = reshape(T.E_phi_imag,   [phiCount, thetaCount]).';

% Reconstruct complex fields
E_theta = E_theta_real + 1i * E_theta_imag;
E_phi   = E_phi_real   + 1i * E_phi_imag;

% Calculations
G_lin = abs(E_theta).^2 + abs(E_phi).^2;
G_dB  = 10 * log10(max(G_lin, 1e-12));

% Circular projection components for Axial Ratio
E_R = (E_theta + 1i * E_phi) / sqrt(2);
E_L = (E_theta - 1i * E_phi) / sqrt(2);
AR_lin = (abs(E_R) + abs(E_L)) ./ max(abs(abs(E_R) - abs(E_L)), 1e-12);
AR_dB  = min(20 * log10(AR_lin), 40); 
AR_dB(G_lin < 1e-4 * max(G_lin(:))) = NaN; % Mask null regions (-40 dB relative)

% Grid construction & 3D Coordinates
[Phi_grid, Theta_grid] = meshgrid(deg2rad(phiVec), deg2rad(thetaVec));
R_plot = max(G_dB - (max(G_dB(:)) - 40), 0);
X = R_plot .* sin(Theta_grid) .* cos(Phi_grid);
Y = R_plot .* sin(Theta_grid) .* sin(Phi_grid);
Z = R_plot .* cos(Theta_grid);

% Rendering Plots
if isFreqDep
    title_suffix = sprintf('%.3f GHz', freqs(f_idx) / 1e9);
else
    title_suffix = 'Frequency Independent';
end

figure('Name', ['Pattern Performance: ' title_suffix], 'Color', 'w', 'Position', [100, 100, 1100, 800]);

% Top-Left: 2D Gain Contour
subplot(2, 2, 1);
contourf(phiVec, thetaVec, G_dB, 25, 'LineColor', 'none');
colorbar; colormap(gca, 'jet'); grid on; clim([max(G_dB(:))-40, max(G_dB(:))]);
title('Total Gain (dB)'); xlabel('\Phi (deg)'); ylabel('\theta (deg)');
set(gca,'YDir','reverse');

% Bottom-Left: 2D Axial Ratio Contour
subplot(2, 2, 3);
contourf(phiVec, thetaVec, AR_dB, 25, 'LineColor', 'none');
colorbar; colormap(gca, 'hot'); grid on; clim([0, 10]);
title('Axial Ratio (dB)'); xlabel('\Phi (deg)'); ylabel('\theta (deg)');
set(gca,'YDir','reverse');

% Top-Right: 3D Spherical Gain
subplot(2, 2, 2);
surf(X, Y, Z, G_dB, 'EdgeColor', 'none'); % 'FaceLighting', 'gouraud'
axis equal; grid on; view(135, 30); % camlight; lighting gouraud; 
colorbar; colormap(gca, 'jet'); clim([max(G_dB(:))-40, max(G_dB(:))]);
title('3D Spherical Gain (dB)'); xlabel('X'); ylabel('Y'); zlabel('Z');

% Bottom-Right: 3D Axial Ratio mapped on Gain
subplot(2, 2, 4);
surf(X, Y, Z, AR_dB, 'EdgeColor', 'none'); % 'FaceLighting', 'gouraud'
axis equal; grid on; view(135, 30); % camlight; lighting gouraud; 
colorbar; colormap(gca, 'hot'); clim([0, 10]);
title('3D Axial Ratio mapped on Gain'); xlabel('X'); ylabel('Y'); zlabel('Z');