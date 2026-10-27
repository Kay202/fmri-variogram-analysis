%  Depth-shell variograms with the partner-validity precondition enforced.
%
%  vg_axis_gamma requires that every x +/- k*e_a stay inside the mask for
%  k <= kmax. Enforcing that is a cross-shaped erosion of radius kmax, so the
%  base set is
%      valid = imerode(eff, vg_cross_elem(KMAX));
%      sel   = valid & depth >= lo & depth < hi;
%  Depth is still measured from the original mask boundary, so the shells keep
%  their meaning; only voxels whose axis partners would fall outside the brain
%  are dropped.
%
%  Same residuals, shells, lags and symmetrisation over +k and -k as the
%  whole-core pass, appended per subject.
%
%  Writes A5_depth_fixed.csv to the output directory set below.

subjects = [3 5 6 8 9 10 11 17 20 21 26 28 31 32 33 34 36 37 38 39 40];
runs     = 1:10;
SHELLS   = [3 4; 4 5; 5 6; 6 8; 8 99];
KMAX     = 3;

here  = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
P = vg_paths('spm', 'fmriprep');

SPMDIR   = P.spm;
RES      = P.residuals;
FMRIPREP = P.fmriprep;
OUT      = fullfile(here, '..', 'data', 'intermediates');
addpath(SPMDIR);
if ~exist(OUT,'dir'), mkdir(OUT); end

fout = fullfile(OUT, 'A5_depth_fixed.csv');
done = [];
if isfile(fout)
    prev = readtable(fout); done = unique(prev.subject)';
    fprintf('already have subjects %s\n', mat2str(done));
end

t0all = tic;
for s = subjects
    if ismember(s, done), fprintf('skipping sub-%02d\n', s); continue; end
    acc = [];
    for r = runs
        if r <= 5, ss = 1; rr = r; else, ss = 2; rr = r - 5; end
        d = fullfile(RES, sprintf('sub%02d', s), sprintf('run%02d', r));
        if isempty(dir(fullfile(d,'Res4D_*.nii')))
            fprintf('  sub-%02d run%02d not packed, skipped\n', s, r); continue;
        end
        t0 = tic;
        [res2, sz, vox, eff] = vg_load(d, FMRIPREP, s, rr, ss);
        depth = bwdist(~eff);
        valid = imerode(eff, vg_cross_elem(KMAX));   % THE FIX

        for sh = 1:size(SHELLS,1)
            raw = eff   & depth >= SHELLS(sh,1) & depth < SHELLS(sh,2);
            sel = valid & depth >= SHELLS(sh,1) & depth < SHELLS(sh,2);
            if nnz(sel) < 500, continue; end
            g = vg_axis_gamma(res2, sz, sel, KMAX);
            for a = 1:3
                for k = 1:KMAX
                    acc = [acc; s, r, ss, SHELLS(sh,1), SHELLS(sh,1)*3, a, k, ...
                           k*vox(a), g(a,k), nnz(sel), nnz(raw)]; %#ok<AGROW>
                end
            end
        end
        clear res2
        fprintf('  sub-%02d run%02d (ses-%02d): %.0f s\n', s, r, ss, toc(t0));
    end
    if ~isempty(acc)
        T = array2table(acc, 'VariableNames', {'subject','run','session', ...
            'shell_lo','depth_mm','axis','k','distance_mm','gamma', ...
            'n_base_valid','n_base_raw'});
        if isfile(fout)
            writetable(T, fout, 'WriteMode','append', 'WriteVariableNames', false);
        else
            writetable(T, fout);
        end
    end
    fprintf('sub-%02d written [%.1f min elapsed]\n', s, toc(t0all)/60);
end
fprintf('\nA5 DONE in %.1f min -> %s\n', toc(t0all)/60, fout);
