clc, clearvars, close all force % Arn_agn ok both

fileName = 'E:\Data\samples\antenna_pattern_sample_HFSS1.ffd';
fileName = 'E:\Data\samples\antenna_pattern_sample_HFSS.ffd';
fileName = 'E:\Data\samples\antenna_pattern_sample_HFSS_2.ffd';

D = ffd_tool(fileName, 'Plot', true);

function D = ffd_tool(filename, varargin)
%FFD_TOOL  Read an HFSS .ffd far-field file; plot Gain & Axial Ratio.
%
%   ffd_tool(file)                read + plot all frequencies
%   ffd_tool(file,'Freq',k)       plot only freq index k
%   ffd_tool(file,'Pin',P)        accepted input power [W] (default 1)
%   ffd_tool(file,'Plot',false)   no plots
%   D = ffd_tool(...)             return struct with all results
%
%   Reads the small header with fgetl/sscanf, then bulk-loads numeric
%   data with readmatrix. Auto-detects Frequency-Independent vs
%   Frequency-Dependent format.
%
%   Phi is normalised to [0,360] (wrap negatives, drop seam duplicate,
%   close sphere at 360) and the data is sorted accordingly.
%
%   Each frequency block is stored as a 6-column TABLE in D.Tables{k}:
%       Theta | Phi | Re_Etheta | Im_Etheta | Re_Ephi | Im_Ephi

% ---------- options ----------
p = inputParser;
p.addParameter('Freq', [],   @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('Pin',  1,    @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('Plot', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});  opt = p.Results;

% ---------- 1) header (3 short reads) ----------
fid = fopen(filename, 'r');
cleaner = onCleanup(@() fclose(fid));
rd = @(s) sscanf(regexprep(s,'[,;]',' '), '%f');     % header-line tokenizer

th = num2cell(rd(fgetl(fid)));   theta = linspace(th{:});   nT = th{3};
ph = num2cell(rd(fgetl(fid)));   phi   = linspace(ph{:});   nP = ph{3};
nF = sscanf(lower(fgetl(fid)), 'frequencies %d');    % [] => independent
isDep        = ~isempty(nF) && nF > 0;
nHeaderLines = 2 + isDep;
nPts         = nT * nP;
clear cleaner

% ---------- 2) frequency values (one regex; readmatrix drops them) ----------
if isDep
    freqs = str2double(regexp(fileread(filename), ...
              '(?<=Frequency\s)[\d.eE+\-]+', 'match', 'ignorecase'));
    assert(numel(freqs) == nF, 'ffd_tool:FreqMismatch', ...
        'Expected %d frequencies, found %d.', nF, numel(freqs));
else
    freqs = NaN;  nF = 1;
end

% ---------- 3) bulk numeric load ("Frequency <v>" rows become all-NaN) ----------
Mfull = readmatrix(filename, 'FileType','text', 'NumHeaderLines',nHeaderLines, ...
                             'ConsecutiveDelimitersRule','join');
MData = Mfull(~any(isnan(Mfull),2), :);              % [nF*nPts x 4]
assert(size(MData,1) == nF*nPts, 'ffd_tool:BadSize', ...
    'Got %d data rows, expected %d (%d freqs × %d pts).', ...
    size(MData,1), nF*nPts, nF, nPts);

% ---------- 4) vectorized split: cell of [nPts x 4] blocks + 3-D tensor ----------
MCell    = mat2cell(MData, repmat(nPts, nF, 1), size(MData,2));     % {nF x 1}
M3D      = permute(reshape(MData, nPts, nF, 4), [1 3 2]);           % [nPts x 4 x nF]
EthetaPF = reshape(complex(M3D(:,1,:), M3D(:,2,:)), nP, nT, nF);    % phi fastest
EphiPF   = reshape(complex(M3D(:,3,:), M3D(:,4,:)), nP, nT, nF);
Et = permute(EthetaPF, [2 1 3]);                                    % [nT x nP x nF]
Ep = permute(EphiPF,   [2 1 3]);

% ---------- 5) PHI CONVENTION: wrap to [0,360], dedupe seam, close sphere ----
phiWrapped = any(phi < 0);
[phi, sortIdx] = unique(mod(phi,360), 'sorted');                    % wrap + dedupe
Et = Et(:, sortIdx, :);
Ep = Ep(:, sortIdx, :);
if phi(1) == 0 && phi(end) < 360                                    % close seam
    phi(end+1)         = 360;
    Et (:, end+1, :)   = Et(:, 1, :);
    Ep (:, end+1, :)   = Ep(:, 1, :);
end
nP   = numel(phi);
nPts = nT * nP;

% ---------- 6) per-block 6-column TABLES (vectorized via cellfun) ----------
%   HFSS row order: (θ1,φ1),(θ1,φ2),...,(θ1,φN),(θ2,φ1),...  -> φ fastest
[PhiCol, ThetaCol] = meshgrid(phi, theta);
ThetaCol = reshape(ThetaCol.', [], 1);
PhiCol   = reshape(PhiCol.',   [], 1);
colNames = {'Theta','Phi','Re_Etheta','Im_Etheta','Re_Ephi','Im_Ephi'};

% Rebuild the sorted/closed per-block [nPts x 4] matrices from Et/Ep
toBlock = @(k) [reshape(real(Et(:,:,k)).', [], 1), reshape(imag(Et(:,:,k)).', [], 1), ...
                reshape(real(Ep(:,:,k)).', [], 1), reshape(imag(Ep(:,:,k)).', [], 1)];
Tables = arrayfun(@(k) array2table([ThetaCol, PhiCol, toBlock(k)], ...
                                   'VariableNames', colNames), ...
                  (1:nF).', 'UniformOutput', false);

% ---------- 7) metrics: vectorized over all frequencies ----------
[PHI, THETA] = meshgrid(phi, theta);
U     = abs(Et).^2 + abs(Ep).^2;
Prad  = reshape(trapz(deg2rad(theta), ...
                trapz(deg2rad(phi), U .* sind(THETA), 2), 1), 1, []);
G_dBi = 10*log10(max(4*pi*U / opt.Pin, eps));
D_dBi = 10*log10(max(4*pi*U ./ reshape(max(Prad,eps),1,1,[]), eps));

ER = (Et - 1j*Ep)/sqrt(2);   EL = (Et + 1j*Ep)/sqrt(2);
AR_dB = 20*log10( (abs(ER)+abs(EL)) ./ max(abs(abs(ER)-abs(EL)), 1e-12) );

% ---------- 8) pack ----------
types = {'independent','dependent'};
D = struct('type',types{1+isDep}, 'nHeaderLines',nHeaderLines, ...
           'phiWrapped',phiWrapped, ...
           'theta_deg',theta, 'phi_deg',phi, 'THETA',THETA, 'PHI',PHI, ...
           'freqs_Hz',freqs, 'Tables',{Tables}, 'MBlocks',{MCell}, ...
           'Etheta',Et, 'Ephi',Ep, ...
           'Gain_dBi',G_dBi, 'Directivity_dBi',D_dBi, ...
           'AxialRatio_dB',AR_dB, 'Prad',Prad, 'filename',filename);

% ---------- 9) plot ----------
if ~opt.Plot, return; end
idx = opt.Freq;  if isempty(idx), idx = 1:nF; end
for i = idx(:).'
    if isnan(freqs(i)), ttl = 'Far Field';
    else,               ttl = sprintf('f = %.4g GHz', freqs(i)/1e9);
    end
    show(THETA, PHI, G_dBi(:,:,i), AR_dB(:,:,i), ttl);
end
end

% =========================================================================
function show(TH, PH, G, AR, ttl)
ARc = min(AR, 40);                                 % clip AR for display
R   = max(G, -30) + 30;                            % gain floor at -30 dBi
sT  = sind(TH);
X = R.*sT.*cosd(PH);  Y = R.*sT.*sind(PH);  Z = R.*cosd(TH);

figure('Name',ttl,'Color','w','Position',[80 80 1100 800]);
ax = gobjects(1,4);

ax(1) = subplot(2,2,1);
contourf(PH,TH,G,20,'LineColor','none'); set(gca,'YDir','reverse');
xlabel('\phi (deg)'); ylabel('\theta (deg)'); colorbar;
title(sprintf('Gain (dBi)  —  peak %.2f', max(G(:))));

ax(2) = subplot(2,2,2);
contourf(PH,TH,ARc,20,'LineColor','none'); set(gca,'YDir','reverse'); hold on;
contour(PH,TH,AR,[3 3],'w-','LineWidth',1.5);
xlabel('\phi (deg)'); ylabel('\theta (deg)'); colorbar;
title(sprintf('Axial Ratio (dB)  —  min %.2f', min(AR(:))));

ax(3) = subplot(2,2,3);
surf(X,Y,Z,G,'EdgeColor','none'); axis equal vis3d off; colorbar;
title('3D Gain (dBi)'); view(135,25); % camlight headlight; lighting gouraud; 

ax(4) = subplot(2,2,4);
surf(X,Y,Z,ARc,'EdgeColor','none'); axis equal vis3d off; colorbar;
title('3D Axial Ratio (dB)'); view(135,25); % camlight headlight; lighting gouraud; 

colormap(ax(1),jet);  colormap(ax(3),jet);
colormap(ax(2),turbo); colormap(ax(4),turbo);
sgtitle(ttl,'FontWeight','bold');
end