function make_stimulus_pair_figure(outPng, smallPath, largePath)
%MAKE_STIMULUS_PAIR_FIGURE  Save a side-by-side stimulus illustration.
%
%   make_stimulus_pair_figure(outPng, smallPath, largePath)
%
% Purpose
%   Convenience utility to display two stimulus images (small vs large) in a
%   clean, publication-ready side-by-side layout and export as a PNG.
%
% Inputs
%   outPng     - Output filename (PNG).
%   smallPath  - Path to the "small / conspecific-sized" stimulus image.
%   largePath  - Path to the "large / predator-sized" stimulus image.
%

% --- read ---
I1 = imread(smallPath);
I2 = imread(largePath);
if size(I1,3)==1, I1 = repmat(I1,[1 1 3]); end
if size(I2,3)==1, I2 = repmat(I2,[1 1 3]); end

% --- match heights for a tidy side-by-side layout ---
h1 = size(I1,1);
h2 = size(I2,1);
targetH = max(h1,h2);
I1s = imresize(I1, [targetH NaN]);
I2s = imresize(I2, [targetH NaN]);

% --- plot ---
f = figure('Color','w','Units','pixels','Position',[100 100 1600 800]);
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

nexttile;
imshow(I1s); axis image off;
title('\bf(A) Conspecific-sized (small)','Interpreter','tex');

nexttile;
imshow(I2s); axis image off;
title('\bf(B) Predator-sized (large)','Interpreter','tex');

% --- save & close ---
[outDir,~,~] = fileparts(outPng);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end
exportgraphics(f, outPng, 'Resolution', 300);
close(f);

end
