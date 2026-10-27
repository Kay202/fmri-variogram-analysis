function Variogram(subjects, runs, outdir)
% Empirical variogram over the Model B residual store, one row per
% subject, run and lag.
%
% The voxel set is a box-eroded core, so an identical set of voxels
% contributes at every lag. Lag vectors are taken over the full half-space by
% explicit indexing, and each lag is symmetrised over +h and -h. Bin
% composition is recorded alongside the estimates.
%
% The per-voxel contribution is computed once per lag and then averaged both
% over the whole core and within depth shells, so the whole-brain curve and
% the depth-resolved curves come out of one pass.
%
% Binning is left to R.

if nargin < 1 || isempty(subjects)
    subjects = [3 5 6 8 9 10 11 17 20 21 26 28 31 32 33 34 36 37 38 39 40];
end
if nargin < 2 || isempty(runs), runs = 1:10; end
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
P = vg_paths('spm', 'fmriprep');

if nargin < 3 || isempty(outdir), outdir = P.out; end

SPMDIR   = P.spm;
RES      = P.residuals;
FMRIPREP = P.fmriprep;
S      = 4;                          % box erosion radius = max lag
SHELLS = [3 4; 4 5; 5 6; 6 8; 8 99]; % depth from boundary, voxels

addpath(SPMDIR);
if ~exist(outdir,'dir'), mkdir(outdir); end

L = halfspace_lags(S);
fprintf('half-space lags: %d\n', size(L,1));

fout = fullfile(outdir, 'rebin_lagwise.csv');
done = [];
if isfile(fout)
    prev = readtable(fout); done = unique(prev.subject)';
    fprintf('already have subjects %s\n', mat2str(done));
end

t0all = tic;
for s = subjects
    if ismember(s, done)
        fprintf('skipping sub-%02d\n', s); continue;
    end
    acc = [];
    for r = runs
        if r <= 5, ss = 1; rr = r; else, ss = 2; rr = r - 5; end
        d = fullfile(RES, sprintf('sub%02d', s), sprintf('run%02d', r));
        if isempty(dir(fullfile(d,'Res4D_*.nii')))
            fprintf('  sub-%02d run%02d not packed, skipped\n', s, r); continue;
        end
        t0 = tic;
        [res2, sz, vox, eff] = vg_load(d, FMRIPREP, s, rr, ss);

        core = imerode(eff, true(2*S+1, 2*S+1, 2*S+1));
        ci = find(core);
        [i1,i2,i3] = ind2sub(sz, ci);
        A = res2(:, ci);
        T = size(res2,1);

        depth = bwdist(~eff);
        dv = depth(ci);
        shell = zeros(numel(ci),1);
        for k = 1:size(SHELLS,1)
            shell(dv >= SHELLS(k,1) & dv < SHELLS(k,2)) = k;
        end

        for q = 1:size(L,1)
            h = L(q,:);
            pv = zeros(numel(ci),1);
            for sgn = [1 -1]
                o = sgn*h;
                si = sub2ind(sz, i1+o(1), i2+o(2), i3+o(3));
                dd = A - res2(:, si);
                pv = pv + 0.5 * (sum(double(dd).^2, 1)' / T);
            end
            pv = pv / 2;                       % symmetrised over +h and -h

            g = nan(1, size(SHELLS,1));
            for k = 1:size(SHELLS,1)
                m = shell == k;
                if nnz(m) >= 100, g(k) = mean(pv(m)); end
            end
            acc = [acc; s, r, ss, h, sqrt(sum((h.*vox).^2)), ...
                   mean(pv), numel(ci), g]; %#ok<AGROW>
        end
        clear res2 A
        fprintf('  sub-%02d run%02d (ses-%02d): %.0f s\n', s, r, ss, toc(t0));
    end
    if ~isempty(acc)
        vn = [{'subject','run','session','dx','dy','dz','distance_mm', ...
               'gamma','n_vox'}, arrayfun(@(k) sprintf('gamma_shell%d',k), ...
               1:size(SHELLS,1), 'uni', 0)];
        Tb = array2table(acc, 'VariableNames', vn);
        if isfile(fout)
            writetable(Tb, fout, 'WriteMode','append', 'WriteVariableNames', false);
        else
            writetable(Tb, fout);
        end
    end
    fprintf('sub-%02d written [%.1f min]\n', s, toc(t0all)/60);
end
fprintf('\nREBIN DONE in %.1f min\n', toc(t0all)/60);
end

function L = halfspace_lags(m)
L = [];
for dx = -m:m, for dy = -m:m, for dz = -m:m
    if dx==0 && dy==0 && dz==0, continue; end
    if dx>0 || (dx==0 && dy>0) || (dx==0 && dy==0 && dz>0)
        L = [L; dx dy dz]; %#ok<AGROW>
    end
end, end, end
end
