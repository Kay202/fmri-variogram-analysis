function VariogramDepth(subjects, runs, outdir)
% Depth-shell and fixed-set variograms over the model B residual store.
%
% Consumes Variogram/residuals_modelB, which holds all 21 subjects across both
% sessions as one float32 4D file per run. Run index follows the authors'
% convention: run01-run05 are ses-01 runs 1-5, run06-run10 are ses-02 runs 1-5,
% and the session matters here because the brain mask lives under the matching
% ses-NN folder.
%
% Writes one row per subject/run/shell/axis/lag so the session comparison (M4)
% and the depth gradient (§6.8) can both be read out in R without recomputing.
%
% Appends per subject, so an interruption resumes rather than discarding work.

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

addpath(SPMDIR);
if ~exist(outdir,'dir'), mkdir(outdir); end

fdep = fullfile(outdir, 'modelB_depth.csv');
fcro = fullfile(outdir, 'modelB_cross.csv');
done = [];
if isfile(fdep)
    prev = readtable(fdep);
    done = unique(prev.subject)';
    fprintf('already have subjects %s\n', mat2str(done));
end

t0all = tic;
for s = subjects
    if ismember(s, done)
        fprintf('skipping sub-%02d (already written)\n', s); continue;
    end
    D = []; C = [];
    for r = runs
        if r <= 5, ss = 1; rr = r; else, ss = 2; rr = r - 5; end
        d = fullfile(RES, sprintf('sub%02d', s), sprintf('run%02d', r));
        if isempty(dir(fullfile(d, 'Res4D_*.nii')))
            fprintf('  sub-%02d run%02d: not packed, skipped\n', s, r);
            continue;
        end
        t0 = tic;
        % vg_run uses `run` only for the mask filename, so pass the within-session
        % index rr; the session itself selects the ses-NN folder
        [Td, Tc] = vg_run(d, FMRIPREP, s, rr, 'B', ss);
        Td.run(:) = r; Tc.run(:) = r;
        Td.session(:) = ss; Tc.session(:) = ss;
        D = [D; Td]; C = [C; Tc]; %#ok<AGROW>
        fprintf('  sub-%02d run%02d (ses-%02d): %.0f s\n', s, r, ss, toc(t0));
    end
    if ~isempty(D)
        append_csv(D, fdep);
        append_csv(C, fcro);
    end
    fprintf('sub-%02d written [%.1f min elapsed]\n', s, toc(t0all)/60);
end
fprintf('\nMODEL B VARIOGRAMS DONE in %.1f min\n', toc(t0all)/60);
end

function append_csv(T, path)
if isempty(T), return; end
if isfile(path)
    writetable(T, path, 'WriteMode','append', 'WriteVariableNames', false);
else
    writetable(T, path);
end
end
