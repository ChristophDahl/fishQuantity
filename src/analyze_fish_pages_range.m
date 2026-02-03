function T = analyze_fish_pages_range(inDir, nStart, nEnd, varargin)
% Batch runner that mirrors test_fish_seg_one (no figures) and computes stats.
% Processes pages nStart..nEnd for sides L and R (default both).
% Writes CSV and overlays beside the input folder.
%
% Extra outputs (per image):
%   AreaTot, AreaMean, AreaStd, PerimTot,
%   HullArea, HullPerim, Density,              % coverage over convex hull
%   CentVar, CentStdX, CentStdY,               % centroid spread across objects
%   NNmean, NNmin, NNmax, PairMean,            % neighbor distances
%   Coverage,                                  % AreaTot / (W*H)
%   MeanEccentricity, MeanSolidity, MeanCompactness  % shape stats
%
% NaN handling:
%   For images with Nobj < 2, NN* and PairMean are undefined. By default
%   we set them to 0 to avoid NaNs in CSV. You can disable with 'FillNNWithZero',false.

% ---------- parameters (IDENTICAL defaults to your tester) ----------
ip = inputParser;
addParameter(ip,'Sides',["L","R"]);
addParameter(ip,'Polarity',"both");      % "both"|"dark"|"bright"
addParameter(ip,'Sensitivity',0.45);     % adaptive threshold sensitivity
addParameter(ip,'Sigma',0.7);            % Gaussian blur sigma
addParameter(ip,'OpenR',1);              % imopen disk radius
addParameter(ip,'CloseR',0);             % imclose disk radius
addParameter(ip,'MinAreaAbs',30);        % absolute minimum area (px)
addParameter(ip,'BandLow',0.40);         % lower keep band as frac*Ahat
addParameter(ip,'BandHigh',2.50);        % upper keep band as frac*Ahat
addParameter(ip,'BigFrac',0.60);         % "big" if Area >= BigFrac*Ahat
addParameter(ip,'DilateFrac',0.6);       % dilation radius ~ frac*sqrt(Ahat)
addParameter(ip,'MaxLinkFactor',1.2);    % link radius ~ factor*sqrt(Ahat)
addParameter(ip,'FillNNWithZero',true);  % replace NN NaNs with 0 when Nobj<2
parse(ip,varargin{:});
P = ip.Results;

inDir = string(inDir);
outCSV = fullfile(inDir, sprintf("run_%04d_%04d_results.csv", nStart, nEnd));
overlayDir = fullfile(inDir, sprintf("run_%04d_%04d_overlays", nStart, nEnd));
if ~exist(overlayDir,'dir'), mkdir(overlayDir); end

% ---------- build file list ----------
files = strings(0,1);
for side = string(P.Sides)
    for n = nStart:nEnd
        f = fullfile(inDir, sprintf("page_%04d_%s.png", n, side));
        if exist(f,'file') == 2
            files(end+1,1) = f; %#ok<AGROW>
        end
    end
end
if isempty(files)
    error('No matching files in %s for %04d..%04d (%s).', inDir, nStart, nEnd, strjoin(string(P.Sides),','));
end

% ---------- pre-alloc ----------
m = numel(files);
Fn = strings(m,1); Page=zeros(m,1); Side=strings(m,1);
W=zeros(m,1); H=zeros(m,1); PolUsed=strings(m,1); Nobj=zeros(m,1);

AreaTot = NaN(m,1); AreaMean = NaN(m,1); AreaStd = NaN(m,1); PerimTot = NaN(m,1);
HullArea = NaN(m,1); HullPerim = NaN(m,1); Density = NaN(m,1);
CentVar = NaN(m,1); CentStdX = NaN(m,1); CentStdY = NaN(m,1);
NNmean = NaN(m,1); NNmin = NaN(m,1); NNmax = NaN(m,1); PairMean = NaN(m,1);
Coverage = NaN(m,1); MeanEccentricity = NaN(m,1); MeanSolidity = NaN(m,1); MeanCompactness = NaN(m,1);

% ---------- main loop ----------
for k = 1:m
    fpath = files(k);
    [~,base,ext] = fileparts(fpath); base = string(base); ext = string(ext);
    Fn(k) = base + ext;

    tok = regexp(char(base), '^page_(\d{4})_([LR])$','tokens','once');
    if ~isempty(tok)
        Page(k) = str2double(tok{1});
        Side(k) = string(tok{2});
    end

    I = imread(fpath);
    if size(I,3)==3, G = rgb2gray(I); else, G = I; end
    G  = im2double(G);
    Gs = imgaussfilt(G, P.Sigma);
    [h,w] = size(G); H(k)=h; W(k)=w;

    if P.Polarity=="both", pols = ["dark","bright"]; else, pols = string(P.Polarity); end

    % evaluate requested polarities; pick best by Nk
    best = struct('Nk',-inf,'BWk',[],'BWfull',[],'pol',"");
    for pol = pols
        % threshold: global vs adaptive -> pick cleaner (lower edge sum)
        BWg = imbinarize(Gs); if pol=="dark", BWg = ~BWg; end
        try
            BWa = imbinarize(Gs,'adaptive','ForegroundPolarity',char(pol),'Sensitivity',P.Sensitivity);
        catch
            Tthr = adaptthresh(Gs, P.Sensitivity, 'ForegroundPolarity', char(pol));
            BWa  = imbinarize(Gs, Tthr);
        end
        s1 = sum(edge(BWg,'sobel'),'all'); s2 = sum(edge(BWa,'sobel'),'all');
        BW = BWg; if s2 < s1, BW = BWa; end

        % light morphology (as in tester)
        if P.OpenR>0,  BW = imopen(BW, strel('disk',P.OpenR)); end
        if P.CloseR>0, BW = imclose(BW, strel('disk',P.CloseR)); end
        BW = imfill(BW,'holes');
        BW = bwareaopen(BW, P.MinAreaAbs);

        % Ahat_pre (pre-merge)
        CC0 = bwconncomp(BW);
        S0  = regionprops(CC0,'Area');
        a0  = vertcat(S0.Area);
        if isempty(a0), Ahat_pre = NaN; else, Ahat_pre = typical_area_from_hist(a0); end

        % robust merge: small -> big (dilate then nearest)
        BWm = merge_small_to_big(BW, Ahat_pre, P.MinAreaAbs, P.BigFrac, P.DilateFrac, P.MaxLinkFactor);

        % recompute + Ahat_post (post-merge)
        CC = bwconncomp(BWm);
        S  = regionprops(CC,'Area','Perimeter','Centroid','Eccentricity','Solidity');
        a  = vertcat(S.Area);

        if isempty(a)
            BWk = false(size(BWm));
            Nk  = 0;
        else
            Ahat_post = typical_area_from_hist(a);
            Ause = Ahat_post; if isnan(Ause), Ause = Ahat_pre; end

            % keep band based on Ause (post-merge preferred)
            minA = max(P.MinAreaAbs, P.BandLow*Ause);
            maxA = P.BandHigh*Ause;
            keep = (a>=minA) & (a<=maxA);
            if ~any(keep) && ~isnan(Ause)
                keep = (a>=0.30*Ause) & (a<=3.00*Ause);
            end

            L   = labelmatrix(CC);
            BWk = ismember(L, find(keep));
            Nk  = nnz(keep);
        end

        if Nk > best.Nk
            best = struct('Nk',Nk,'BWk',BWk,'BWfull',BWm,'pol',string(pol));
        end
    end

    % store + overlay
    BWk = best.BWk; BWfull = best.BWfull; Nobj(k) = best.Nk; PolUsed(k) = best.pol;

    try
        rgb = I; if size(rgb,3)==1, rgb = repmat(I,[1 1 3]); end
        perK = bwperim(BWk);              [ry,rx] = find(perK); idx = sub2ind([h w],ry,rx);
        rgb(idx) = 0; rgb(idx+h*w) = 255; rgb(idx+2*h*w) = 0;  % green
        perR = bwperim(BWfull & ~BWk);    [ry,rx] = find(perR); idx = sub2ind([h w],ry,rx);
        rgb(idx) = 255; rgb(idx+h*w) = 0; rgb(idx+2*h*w) = 0;  % red
        imwrite(rgb, char(fullfile(overlayDir, base + ext)));
    catch
        % ignore overlay errors
    end

    % ---------- statistics on kept objects ----------
    if Nobj(k) > 0
        CCk = bwconncomp(BWk);
        Sk  = regionprops(CCk,'Area','Perimeter','Centroid','Eccentricity','Solidity');
        ak = vertcat(Sk.Area);
        pk = vertcat(Sk.Perimeter);
        Ck = vertcat(Sk.Centroid);

        % areas / perimeters
        AreaTot(k) = sum(ak);
        AreaMean(k)= mean(ak);
        AreaStd(k) = std(ak);
        PerimTot(k)= sum(pk);

        % convex hull & density (using all kept pixels)
        [yy,xx] = find(BWk);
        if numel(xx) >= 3
            K = convhull(xx,yy);
            HullArea(k)  = polyarea(xx(K), yy(K));
            HullPerim(k) = sum(hypot(diff(xx(K)), diff(yy(K))));
            Density(k)   = AreaTot(k) / max(HullArea(k), eps);
        else
            HullArea(k) = NaN; HullPerim(k) = NaN; Density(k) = NaN;
        end

        % page coverage
        Coverage(k) = AreaTot(k) / (W(k)*H(k));

        % shape descriptors
        ecc = vertcat(Sk.Eccentricity);
        sol = vertcat(Sk.Solidity);
        MeanEccentricity(k) = mean(ecc);
        MeanSolidity(k)     = mean(sol);
        MeanCompactness(k)  = mean( 4*pi*ak ./ max(pk.^2, eps) );

        % centroid spread & neighbor distances
        if size(Ck,1) >= 2
            Ccov = cov(Ck,1);                   % population covariance
            CentVar(k) = trace(Ccov);
            CentStdX(k) = sqrt(Ccov(1,1));
            CentStdY(k) = sqrt(Ccov(2,2));

            dvec = pdist(Ck);                    % pairwise distances (vector form)
            PairMean(k) = mean(dvec);

            Dm = squareform(dvec);
            Dm(Dm==0) = NaN;
            nn = min(Dm,[],2,'omitnan');
            NNmean(k) = mean(nn,'omitnan');
            NNmin(k)  = min(nn,[],'omitnan');
            NNmax(k)  = max(nn,[],'omitnan');
        else
            % single object: spread is zero; NN undefined -> fill per setting
            CentVar(k)  = 0;
            CentStdX(k) = 0;
            CentStdY(k) = 0;
            if P.FillNNWithZero
                NNmean(k) = 0; NNmin(k) = 0; NNmax(k) = 0; PairMean(k) = 0;
            else
                % leave as NaN
            end
        end
    else
        % no objects: set everything to 0 for stability (optional)
        if P.FillNNWithZero
            CentVar(k)=0; CentStdX(k)=0; CentStdY(k)=0;
            NNmean(k)=0; NNmin(k)=0; NNmax(k)=0; PairMean(k)=0;
            Coverage(k)=0; Density(k)=0;
            AreaTot(k)=0; AreaMean(k)=0; AreaStd(k)=0; PerimTot(k)=0;
            MeanEccentricity(k)=0; MeanSolidity(k)=0; MeanCompactness(k)=0;
            HullArea(k)=0; HullPerim(k)=0;
        end
    end
end

T = table(Fn,Page,Side,W,H,PolUsed,Nobj, ...
          AreaTot,AreaMean,AreaStd,PerimTot, ...
          HullArea,HullPerim,Density,Coverage, ...
          CentVar,CentStdX,CentStdY, ...
          NNmean,NNmin,NNmax,PairMean, ...
          MeanEccentricity,MeanSolidity,MeanCompactness);

writetable(T, char(outCSV));
fprintf('Saved results to: %s\nOverlays in: %s\n', outCSV, overlayDir);
end

% ================= helpers =================
function BWm = merge_small_to_big(BW, Ahat, MinAreaAbs, BigFrac, DilateFrac, MaxLinkFactor)
% Merge each "small" component into a nearby "big" one (touch via dilation,
% then nearest-centroid fallback). Safe no-op if Ahat is invalid.
if isempty(Ahat) || isnan(Ahat) || Ahat<=0, BWm = BW; return; end
CC = bwconncomp(BW);
S  = regionprops(CC,'Area','Centroid');
if isempty(S), BWm = BW; return; end
a = vertcat(S.Area);
C = vertcat(S.Centroid);
bigCut = max(MinAreaAbs, BigFrac*Ahat);
isBig  = a >= bigCut;
isSmall= ~isBig;
if ~any(isBig)                                   % anchor if all small
    [~,imax] = max(a); isBig(imax)=true; isSmall(imax)=false;
end
L   = labelmatrix(CC);
BWm = false(size(BW));
for j = find(isBig).'
    BWm = BWm | (L==j);
end
R  = max(1, round(DilateFrac * sqrt(Ahat)));     % touch-based merge
SE = strel('disk', R);
BWd = imdilate(BWm, SE);
smallIdx = find(isSmall);
assigned = false(size(smallIdx));
for s = 1:numel(smallIdx)
    Msm = (L==smallIdx(s));
    if any(BWd & Msm, 'all')
        BWm = BWm | Msm; assigned(s) = true;
    end
end
left = smallIdx(~assigned);                      % nearest-centroid fallback
if ~isempty(left)
    Cbig = C(isBig,:); Csm = C(left,:);
    D = pdist2(Csm, Cbig);
    [dmin, ~] = min(D,[],2);
    MaxLink = MaxLinkFactor * sqrt(Ahat);
    for k = 1:numel(left)
        if isfinite(dmin(k)) && dmin(k) <= MaxLink
            BWm = BWm | (L==left(k));
        end
    end
end
BWm = imclose(BWm, strel('disk',1));            % cleanup
BWm = imfill(BWm,'holes');
BWm = bwareaopen(BWm, MinAreaAbs);
end

function A = typical_area_from_hist(a)
% Robust “typical area”: mean of mode bin; fallback to median.
if isempty(a), A = NaN; return; end
nb = min(50, max(15, round(sqrt(numel(a)))));
[cnt,edges] = histcounts(a, nb);
if any(cnt)
    [~,imx] = max(cnt);
    A = mean(edges([imx imx+1]));
else
    A = median(a);
end
end
