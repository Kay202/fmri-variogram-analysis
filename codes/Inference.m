% The cluster-extent inference experiment, in four stages:
%
%   1. reference fidelity   does the analytic reference synthesise faithfully
%   2. error rates          familywise error for each assumed model
%   3. thresholds           what the unidentified models imply for c*
%   4. power                detection at the calibrated thresholds
%
% Run Analysis.R first: this reads S25_reference_params.csv, the S27_ref_* fits
% and the S34_* sweep.
%
% Each stage writes one CSV and is skipped if that CSV already exists, so an
% interrupted run continues where it stopped. Delete a CSV to force a re-run.
%
% Needs the Image Processing Toolbox for bwconncomp.

restoredefaultpath; rehash toolboxcache;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'helpers'));
OUTV = fullfile(here, '..', 'data', 'intermediates');
if ~exist(OUTV, 'dir'), mkdir(OUTV); end

stages = {
    'S26_reference_fidelity.csv', @stage_S26, 'does the analytic reference synthesise faithfully'
    'S32_e1e10_fwer_final.csv', @stage_S32, 'the six-arm error-rate table'
    'S35_sweep_thresholds.csv', @stage_S35, 'what the unidentified models imply for thresholds'
    'S29_e1e10_power_ref.csv', @stage_S29, 'power at the calibrated thresholds'
};

fprintf('inference experiment: %d stages\n\n', size(stages, 1));
t_all = tic;
for i = 1:size(stages, 1)
    f = fullfile(OUTV, stages{i, 1});
    if exist(f, 'file') == 2
        fprintf('[have] %-28s %s\n', stages{i, 1}, stages{i, 3});
        continue
    end
    fprintf('[run ] %-28s %s\n', stages{i, 1}, stages{i, 3});
    t0 = tic;
    stages{i, 2}(OUTV);
    if exist(f, 'file') ~= 2
        error('%s finished without writing %s', func2str(stages{i, 2}), stages{i, 1});
    end
    fprintf('       %.1f min\n', toc(t0)/60);
end
fprintf('\ninference experiment complete, %.1f min total\n', toc(t_all)/60);

% 1. does the analytic reference synthesise faithfully
function stage_S26(OUTV)
%% Stage 1. Fidelity of the analytic reference
%
%  Synthesises a field from the reference covariance and compares the
%  variogram it realises with the reference's own, so the rest of the
%  experiment is only run if the generator is faithful. The reference has a
%  spectral density in closed form, evaluated directly at every mode:
%
%      f(k)  proportional to  sum_m s_m tau_m
%                             ( exp(-(a_m-k)^2 tau_m^2/2)
%                             - exp(-(a_m+k)^2 tau_m^2/2) ) / (k a_m)
%
%  Writes S26_reference_fidelity.csv to the intermediates directory.


vox = 3; gridsz = [64 64 64];

kmax = 4; V = [];
for dx = 0:kmax, for dy = 0:kmax, for dz = 0:kmax
    if dx == 0 && dy == 0 && dz == 0, continue; end
    V = [V; dx dy dz dx^2+dy^2+dz^2]; %#ok<AGROW>
end, end, end
k2u = unique(V(:,4)); k2u = k2u(k2u <= 48);

P = readtable(fullfile(OUTV, 'S25_reference_params.csv'));
Kmag = mode_freqs(gridsz, vox);
fprintf('grid %d mm, fundamental %.4f rad/mm, Nyquist corner %.4f rad/mm\n', ...
        gridsz(1)*vox, min(Kmag(Kmag>0)), max(Kmag(:)));

rows = [];
for w = {'matched','hole'}
    which = w{1};
    fprintf('\n===== %s reference =====\n', which);
    fprintf('  %-6s %10s %12s %12s %14s\n', 'shell', 'modes', 'max err', 'rms err', 'as %% of range');
    for k = [3 4 5]
        t = P(strcmp(P.reference, which) & P.shell == k, :);
        s = t.sill; a = t.a; tau = t.tau; nug = t.nugget(1);

        % the variogram the reference specifies
        g_exact = zeros(numel(k2u),1);
        for j = 1:numel(k2u)
            r = vox*sqrt(k2u(j));
            g_exact(j) = nug + sum(s .* (1 - om3(a*r).*exp(-r^2./(2*tau.^2))));
        end

        % the spectrum the synthesiser will use, and the field it implies
        S = spec_ref(Kmag, s, a, tau, sum(s));
        g_gen = gen_gamma(S.^2, nug, V, k2u, gridsz);
        err = g_gen - g_exact;
        rng_ = max(g_exact) - min(g_exact);
        fprintf('  %-6d %10d %12.4f %12.4f %13.2f%%\n', k, nnz(S > 0), ...
                max(abs(err)), sqrt(mean(err.^2)), 100*max(abs(err))/rng_);
        rows = [rows; k, strcmp(which,'hole'), nnz(S>0), max(abs(err)), ...
                sqrt(mean(err.^2)), rng_]; %#ok<AGROW>
    end
end

R = array2table(rows, 'VariableNames', ...
    {'shell','is_hole','modes_used','max_err','rms_err','curve_range'});
writetable(R, fullfile(OUTV, 'S26_reference_fidelity.csv'));

fprintf('\n================ verdict ================\n');
fprintf('  the four candidate models differ from each other by at most 0.023\n');
fprintf('  the shell-binned generator missed its own model by 0.134 to 0.429\n');
fprintf('  the analytic reference misses its own model by %.4f at worst (%.2f%% of range)\n', ...
        max(R.max_err), 100*max(R.max_err./R.curve_range));
if max(R.max_err) < 0.023
    fprintf('  Faithful to inside the differences being measured. Proceed.\n');
else
    fprintf('  NOT faithful enough. Do not run the experiment on this.\n');
end

%% ========================================================================
end

% 2. the six-arm error-rate table
function stage_S32(OUTV)
%% Stage 2. Familywise error rate, six arms
%
%  Draws null fields from the reference covariance, calibrates a critical
%  cluster size for each assumed model, and scores every threshold against
%  that one truth.
%
%  The six arms are the four assumed models, the admissible estimate in two
%  forms, and the reference itself as a labelled check that the machinery
%  returns the nominal rate when the assumed model is correct. The two
%  admissible arms differ only in how fields are drawn from the estimate: as
%  normally synthesised, and with each spectral shell damped,
%
%      C(r) = p0*1{r=0} + sum_j p_j Omega_3(t_j r) exp(-r^2 / (2 tau^2))
%
%  which keeps every term a valid covariance in R^3 and spreads each shell
%  over enough lattice modes to be representable. tau = 80 mm.
%
%  Writes S32_e1e10_fwer_final.csv to the intermediates directory.


nsim = 1000; gridsz = [64 64 64]; vox = 3; zthr = [3.09 2.30];
TAU_DAMP = 80;


REF = readtable(fullfile(OUTV, 'S25_reference_params.csv'));
Kmag = mode_freqs(gridsz, vox);
names = {'Gaussian','AFNImixed','Exponential','SB','SBdamped','Reference'};
rows = [];

for w = {'matched','hole'}
    which = w{1};
    C = readtable(fullfile(OUTV, ['S27_ref_' which '_sbcoefs.csv']));
    P = readtable(fullfile(OUTV, ['S27_ref_' which '_params.csv']));
    E = readtable(fullfile(OUTV, ['S27_ref_' which '_exp.csv']));
    fprintf('\n############ %s reference ############\n', which);

    for si = 1:height(P)
        k = P.shell(si); dmm = P.depth_mm(si);
        ei = find(E.shell == k, 1);
        t  = REF(strcmp(REF.reference, which) & REF.shell == k, :);
        fprintf('\n===== depth shell %d (%d mm) =====\n', k, dmm);

        c = C(C.shell == k, :); p_sb = c.coef; nodes = c.node(2:end);
        act = p_sb(2:end) > 1e-12; tk = nodes(act);
        pk = p_sb(2:end); pk = pk(act);

        A = struct(); nug = struct();
        A.Reference   = spec_ref(Kmag, t.sill, t.a, t.tau, sum(t.sill));
        nug.Reference = t.nugget(1);
        A.truth = A.Reference;  nug_truth = t.nugget(1);
        A.SB          = vg_spectral_field(p_sb, nodes, gridsz, vox);
        nug.SB        = p_sb(1);
        A.SBdamped    = spec_ref(Kmag, pk, tk, repmat(TAU_DAMP, numel(tk), 1), sum(pk));
        nug.SBdamped  = p_sb(1);
        A.Gaussian    = spec_gauss(Kmag, P.g_range(si), P.g_sill(si));
        nug.Gaussian  = P.g_nugget(si);
        A.AFNImixed   = spec_mixed(Kmag, P.m_b(si), P.m_c(si), P.m_p1(si), P.m_p2(si));
        nug.AFNImixed = 0;
        A.Exponential = spec_exp(Kmag, E.e_range(ei), E.e_sill(ei) - E.e_nugget(ei));
        nug.Exponential = E.e_nugget(ei);

        for ti = 1:numel(zthr)
            z = zthr(ti);
            fprintf('  -- cluster-forming threshold z > %.2f --\n', z);
            cstar = zeros(1, numel(names)); gam = zeros(1, numel(names));
            for m = 1:numel(names)
                mx = max_clusters(A.(names{m}), nug.(names{m}), gridsz, nsim, z, ...
                                  4000 + 10*ti + m);
                [cstar(m), gam(m)] = crit_randomised(mx, 0.05);
                fprintf(['    %-12s assumed: c* = %5.0f voxels, tie prob %.2f ' ...
                         '(median max %5.0f, %d distinct)\n'], ...
                        names{m}, cstar(m), gam(m), median(mx), numel(unique(mx)));
            end
            mxT = max_clusters(A.truth, nug_truth, gridsz, nsim, z, 7000 + ti);
            fprintf('    truth       : median max cluster %5.0f\n', median(mxT));
            rng(9000 + ti);
            for m = 1:numel(names)
                fwer_ge = mean(mxT >= cstar(m));
                u = rand(size(mxT));
                fwer = mean( (mxT > cstar(m)) | (mxT == cstar(m) & u < gam(m)) );
                se = sqrt(fwer * (1 - fwer) / nsim);
                fprintf(['    --> assuming %-12s : FWER = %.3f  (+/- %.3f)   ' ...
                         '[inclusive rule %.3f]\n'], names{m}, fwer, 1.96*se, fwer_ge);
                rows = [rows; strcmp(which,'hole'), k, dmm, z, m, cstar(m), ...
                        gam(m), fwer, 1.96*se, fwer_ge]; %#ok<AGROW>
            end
        end
    end
end

T = array2table(rows, 'VariableNames', {'is_hole','shell','depth_mm','zthr', ...
    'model','c_star','gamma','FWER','ci95','FWER_inclusive'});
T.model_name = repmat(names(:), height(T)/numel(names), 1);
T.reference = repmat({'matched'}, height(T), 1);
T.reference(T.is_hole == 1) = {'hole'};
writetable(T, fullfile(OUTV, 'S32_e1e10_fwer_final.csv'));

fprintf('\n================ summary, nominal FWER = 0.05 ================\n');
for w = {'matched','hole'}
    fprintf('  --- %s reference ---\n', w{1});
    for m = 1:numel(names)
        sel = T.model == m & strcmp(T.reference, w{1});
        fprintf('    %-12s c* %5.0f to %5.0f   FWER %.3f to %.3f  (mean %.3f)\n', ...
                names{m}, min(T.c_star(sel)), max(T.c_star(sel)), ...
                min(T.FWER(sel)), max(T.FWER(sel)), mean(T.FWER(sel)));
    end
end
fprintf('\nReference is the machinery check. SB against SBdamped separates the\n');
fprintf('estimator from the way fields are drawn from it.\n');
fprintf('\nwrote S32_e1e10_fwer_final.csv (%d simulations per cell)\n', nsim);

%% ===== helpers =============================================================
end

% 3. what the unidentified models imply for thresholds
function stage_S35(OUTV)
%% Stage 3. Cluster thresholds across the unidentified models
%
%  Takes the models the identifiability sweep admitted, synthesises each one
%  analytically through the closed-form density, and records the critical
%  cluster size each would hand a practitioner. Every threshold is scored
%  against one fixed truth, the matched reference, so the spread in realised
%  error rate is attributable to a choice the data cannot make.
%
%  The admission tolerance is recorded per row rather than fixed, so the
%  result can be read at whichever criterion is preferred.
%
%  Writes S35_sweep_thresholds.csv to the intermediates directory.


nsim = 1000; gridsz = [64 64 64]; vox = 3; zthr = [3.09 2.30];

P   = readtable(fullfile(OUTV, 'S34_sweep_params.csv'));
REF = readtable(fullfile(OUTV, 'S25_reference_params.csv'));
S34 = readtable(fullfile(OUTV, 'S34_sweep.csv'));
Kmag = mode_freqs(gridsz, vox);

rows = [];
for k = [3 4 5]
    floor_k = min(S34.rmse(S34.shell == k));
    tol_k   = S34.tol(find(S34.shell == k, 1));
    keep = unique(P.total(P.shell == k & P.rmse <= tol_k));
    t = REF(strcmp(REF.reference, 'matched') & REF.shell == k, :);
    fprintf('\n===== shell %d : %d models within the curve''s own SE =====\n', ...
            k, numel(keep));

    Strue = spec_ref(Kmag, t.sill, t.a, t.tau, sum(t.sill));
    nug_true = t.nugget(1);

    for ti = 1:numel(zthr)
        z = zthr(ti);
        mxT = max_clusters(Strue, nug_true, gridsz, nsim, z, 7000 + ti);
        fprintf('  -- z > %.2f, truth median max cluster %d --\n', z, median(mxT));
        fprintf('  %8s %9s %9s %8s %8s\n', 'total', 'rmse', 'c*', 'FWER', 'within');
        rng(9500 + ti);
        for j = 1:numel(keep)
            q = P(P.shell == k & abs(P.total - keep(j)) < 1e-9, :);
            S = spec_ref(Kmag, q.sill, q.a, q.tau, sum(q.sill));
            mx = max_clusters(S, q.nugget(1), gridsz, nsim, z, 5000 + 100*ti + j);
            [cs, gm] = crit_randomised(mx, 0.05);
            u = rand(size(mxT));
            fwer = mean( (mxT > cs) | (mxT == cs & u < gm) );
            lab = 'SE'; if q.rmse(1) <= 3*floor_k, lab = '3x floor'; end
            if q.rmse(1) <= 1.5*floor_k, lab = '1.5x floor'; end
            fprintf('  %8.3f %9.4f %9.0f %8.3f %8s\n', keep(j), q.rmse(1), cs, fwer, lab);
            rows = [rows; k, z, keep(j), q.rmse(1), q.nugget(1), cs, gm, fwer, ...
                    floor_k, tol_k]; %#ok<AGROW>
        end
    end
end

T = array2table(rows, 'VariableNames', {'shell','zthr','total','rmse','nugget', ...
    'c_star','gamma','FWER','floor_rmse','tol'});
writetable(T, fullfile(OUTV, 'S35_sweep_thresholds.csv'));

fprintf('\n================ the bound ================\n');
fprintf('Every model below reproduces the observed curve to within the stated\n');
fprintf('tolerance, so the data cannot choose among them. Thresholds are scored\n');
fprintf('against one fixed truth, so the spread is what the ambiguity costs.\n\n');
crits = {1.5, '1.5x the best achievable fit'; 3, '3x the best achievable fit'; ...
         Inf, 'the curve''s own standard error'};
fprintf('  %-32s %6s %6s %14s %16s\n', 'tolerance', 'shell', 'z', 'c* range', 'FWER range');
for c = 1:size(crits,1)
    for k = [3 4 5]
        for z = zthr
            s = T(T.shell == k & abs(T.zthr - z) < 1e-9, :);
            if isfinite(crits{c,1}), s = s(s.rmse <= crits{c,1}*s.floor_rmse, :); end
            if isempty(s), continue; end
            fprintf('  %-32s %6d %6.2f %6.0f to %-4.0f %7.3f to %-6.3f\n', ...
                    crits{c,2}, k, z, min(s.c_star), max(s.c_star), ...
                    min(s.FWER), max(s.FWER));
        end
    end
end
fprintf('\nA nominal rate of 0.05 is the target throughout.\n');
fprintf('wrote S35_sweep_thresholds.csv (%d simulations per cell)\n', nsim);

%% ========================================================================
end

% 4. power at the calibrated thresholds
function stage_S29(OUTV)
%% Stage 4. Power at the calibrated thresholds
%
%  Adds a signal to fields drawn from the reference covariance and measures
%  detection at the threshold each method produced in stage 2, read alongside
%  the error rate that threshold realises.
%
%  Two details that are easy to get wrong. Thresholds are selected by name
%  rather than by model index, because the index order differs between runs.
%  And the tie probability travels with the critical value: the null
%  distribution of the maximum cluster is discrete and the calibration uses a
%  randomised rule, so power must be computed under the same rule.
%
%  Reads  S25_reference_params.csv, S27_ref_<which>_{sbcoefs,params,exp}.csv
%         and the error rates from stage 2, all in the intermediates directory
%  Writes S29_e1e10_power_ref.csv to the intermediates directory.


nsim   = 500;
amps   = [1 2 3 4 5 6];
gridsz = [64 64 64];
vox    = 3;
zthr   = 3.09;
FWHM = 15;   % a 9 mm blob's suprathreshold core is smaller than the thresholds


REF = readtable(fullfile(OUTV, 'S25_reference_params.csv'));
% calibration comes from the six-arm run, so the damped admissible estimate is
% carried through to power as well; S28's five-arm table is superseded
F   = readtable(fullfile(OUTV, 'S32_e1e10_fwer_final.csv'));
F   = F(abs(F.zthr - zthr) < 1e-6, :);
Kmag = mode_freqs(gridsz, vox);

sig  = blob(gridsz, vox, FWHM);
mask = sig > 0.5 * max(sig(:));
fprintf('signal blob: FWHM %g mm, %d voxels above half-max\n', FWHM, nnz(mask));

names = {'Gaussian','AFNImixed','Exponential','SB','SBdamped','Reference'};
nm = numel(names); rows = [];

for w = {'matched','hole'}
    which = w{1};
    P = readtable(fullfile(OUTV, ['S27_ref_' which '_params.csv']));
    fprintf('\n############ %s reference ############\n', which);
    for si = 1:height(P)
        k = P.shell(si); dmm = P.depth_mm(si);
        t = REF(strcmp(REF.reference, which) & REF.shell == k, :);
        S = spec_ref(Kmag, t.sill, t.a, t.tau, sum(t.sill));
        nug = t.nugget(1); N = prod(gridsz);

        cs = zeros(1,nm); gm = zeros(1,nm); fw = zeros(1,nm);
        for m = 1:nm
            r = F(F.shell == k & strcmp(F.model_name, names{m}) & ...
                  strcmp(F.reference, which), :);
            if isempty(r)
                error('no calibration row: shell %d, %s, %s', k, names{m}, which);
            end
            cs(m) = r.c_star(1); gm(m) = r.gamma(1); fw(m) = r.FWER(1);
        end
        fprintf('\n===== shell %d (%d mm), z > %.2f =====\n', k, dmm, zthr);
        for m = 1:nm
            fprintf('  %-12s threshold %6.0f voxels, realised FWER %.3f\n', ...
                    names{m}, cs(m), fw(m));
        end

        for a = amps
            rng(round(1e4*a) + k + 500*strcmp(which,'hole'));
            det = zeros(1,nm);
            for i = 1:nsim
                f = real(ifftn(S .* (randn(gridsz) + 1i*randn(gridsz)))) * N;
                if nug > 0, f = f + randn(gridsz) * sqrt(nug); end
                s0 = std(f(:));
                f = f + a * s0 * sig;
                z = (f - mean(f(:))) / std(f(:));
                cc = bwconncomp(z > zthr, 26);
                best = 0;
                for j = 1:cc.NumObjects
                    idx = cc.PixelIdxList{j};
                    if any(mask(idx)), best = max(best, numel(idx)); end
                end
                det = det + ((best > cs) | (best == cs & rand(1,nm) < gm));
            end
            pw = det / nsim;
            fprintf('    amplitude %.1f SD : %s\n', a, ...
                    strjoin(arrayfun(@(m) sprintf('%s %.3f', names{m}, pw(m)), ...
                                     1:nm, 'UniformOutput', false), ' | '));
            for m = 1:nm
                rows = [rows; strcmp(which,'hole'), k, dmm, zthr, m, a, ...
                        cs(m), fw(m), pw(m)]; %#ok<AGROW>
            end
        end
    end
end

T = array2table(rows, 'VariableNames', {'is_hole','shell','depth_mm','zthr', ...
    'model','amplitude','c_star','FWER','power'});
T.model_name = repmat(names(:), height(T)/nm, 1);
T.reference = repmat({'matched'}, height(T), 1);
T.reference(T.is_hole == 1) = {'hole'};
writetable(T, fullfile(OUTV, 'S29_e1e10_power_ref.csv'));

fprintf('\n================ summary ================\n');
fprintf('Only arms holding nominal error have interpretable power.\n');
for w = {'matched','hole'}
    fprintf('\n  --- %s reference ---\n', w{1});
    for m = 1:nm
        sel = T.model == m & strcmp(T.reference, w{1});
        fprintf('    %-12s realised FWER %.3f, power at amp 3 = %.3f\n', names{m}, ...
                mean(T.FWER(sel)), mean(T.power(sel & T.amplitude == 3)));
    end
end
fprintf('\nwrote S29_e1e10_power_ref.csv (%d simulations per cell)\n', nsim);

%% ------------------------------------------------------------------------
% reproduced verbatim from Diagnostics/matlab/e1e10_power.m, where it is a
% local function and so not visible outside that file
end

% helpers used by one stage only

function S = spec_gauss(K, a, sill)
P = exp(-(K.^2) * a^2 / 4);  P(1) = 0;
S = sqrt(P * sill / sum(P(:)));
end


function S = spec_mixed(K, b, c, p1, p2)
Pg = exp(-(K.^2) * b^2 / 2);      Pg(1) = 0;
Pe = 1 ./ (1 + (K * c).^2).^2;    Pe(1) = 0;
P = p1 * Pg / sum(Pg(:)) + p2 * Pe / sum(Pe(:));
S = sqrt(P);
end


function S = spec_exp(K, phi, sill)
P = 1 ./ (1 + (K * phi).^2).^2;  P(1) = 0;
S = sqrt(P * sill / sum(P(:)));
end

function b = blob(gridsz, vox, fwhm)
c = (gridsz + 1) / 2;
[X,Y,Z] = ndgrid(1:gridsz(1), 1:gridsz(2), 1:gridsz(3));
r2 = ((X-c(1)).^2 + (Y-c(2)).^2 + (Z-c(3)).^2) * vox^2;
s = fwhm / (2*sqrt(2*log(2)));
b = exp(-r2 / (2*s^2));
end

%% ========================================================================

