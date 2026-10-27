% Why the reference covariance is analytic rather than shell-binned: the shell
% binner's synthesised field misses its own model by more than the differences
% the inference experiment measures.
%
% Three independent blocks. Produces no manuscript table.
%
% Reads and writes the intermediate result tables, in the directory set below.

restoredefaultpath; rehash toolboxcache;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'helpers'));
OUTV = fullfile(here, '..', 'data', 'intermediates');

% Block 1. Fidelity of the shell-binned generator to its own model
%
%  Builds a field from each of four admissible models with the shell-binning
%  synthesiser, and compares the variogram realised against the model's own.
%
%  The synthesiser assigns each node's mass to the Fourier modes in a shell,
%  with edges at the midpoints between active node radii, so a shell of zero
%  thickness cannot be represented on a lattice whose modes are spaced
%  2*pi/L apart. The departure reported here is that effect.


gridsz = [64 64 64]; vox = 3;

% the 31 lattice separations and one representative displacement for each
kmax = 4; V = [];
for dx = 0:kmax, for dy = 0:kmax, for dz = 0:kmax
    if dx == 0 && dy == 0 && dz == 0, continue; end
    V = [V; dx dy dz dx^2+dy^2+dz^2]; %#ok<AGROW>
end, end, end
k2u = unique(V(:,4)); k2u = k2u(k2u <= 48);
fprintf('%d distinct separations from %.2f to %.2f mm\n', ...
        numel(k2u), vox*sqrt(min(k2u)), vox*sqrt(max(k2u)));

models = {'SB','NCS','LP','PTR'};
files  = {'shell_sb_coefs_v2.csv','S20_coefs_NCS.csv','S20_coefs_LP.csv','S20_coefs_PTR.csv'};
T = struct();
for i = 1:numel(models), T.(models{i}) = readtable(fullfile(OUTV, files{i})); end

rows = [];
for k = [3 4 5]
    fprintf('\n===== shell %d =====\n', k);
    fprintf('  %-5s %12s %12s %14s %14s\n', 'model', 'lowest bin', 'modes in it', ...
            'max gamma err', 'rms gamma err');
    for i = 1:numel(models)
        t = T.(models{i}); t = t(t.shell == k, :);
        p = t.coef; nodes = t.node(2:end);
        [S, info] = vg_spectral_field(p, nodes, gridsz, vox);
        P = S.^2;

        % the covariance the generator will actually produce
        C = real(ifftn(P)) * prod(gridsz);

        % the variogram the model specifies
        w = p(2:end); tk = nodes(w > 0); wk = w(w > 0);
        g_exact = zeros(numel(k2u),1); g_gen = zeros(numel(k2u),1);
        for j = 1:numel(k2u)
            r = vox*sqrt(k2u(j));
            g_exact(j) = p(1) + sum(wk .* (1 - sin(tk*r)./(tk*r)));
            % average the generator's covariance over every displacement of this length
            sel = V(V(:,4) == k2u(j), 1:3);
            acc = 0;
            for m = 1:size(sel,1)
                idx = mod(sel(m,:), gridsz) + 1;
                acc = acc + C(idx(1), idx(2), idx(3));
            end
            g_gen(j) = p(1) + (C(1,1,1) - acc/size(sel,1));
        end
        err = g_gen - g_exact;

        active = sort(tk);
        lowbin = (active(1) + active(2))/2;
        nmodes = info.modes_used;
        fprintf('  %-5s %12.4f %12d %14.4f %14.4f\n', models{i}, lowbin, nmodes, ...
                max(abs(err)), sqrt(mean(err.^2)));
        rows = [rows; k, i, lowbin, nmodes, max(abs(err)), sqrt(mean(err.^2)), ...
                info.frac_mass_lost]; %#ok<AGROW>
    end
end

R = array2table(rows, 'VariableNames', ...
    {'shell','model','low_bin_edge','modes_used','max_err','rms_err','frac_mass_lost'});
R.model_name = repmat(models(:), height(R)/numel(models), 1);
writetable(R, fullfile(OUTV, 'S21_generator_fidelity.csv'));

fprintf('\n================ verdict ================\n');
fprintf('  curve ranges are about 1.5; the four models agree with each other to 0.023\n');
fprintf('  worst generator error across all models and shells: %.4f\n', max(R.max_err));
if max(R.max_err) > 0.023
    fprintf('  The generator departs from its own model by MORE than the models\n');
    fprintf('  differ from each other. The S20 spread is manufactured.\n');
else
    fprintf('  The generator is faithful to well inside the model differences.\n');
    fprintf('  The S20 spread is real.\n');
end

% Block 2. Whether the departure is in the model or in the synthesis
%
%  Separates the two by comparing, for the admissible estimate, the variogram
%  of its own curve against the variogram of a field drawn from it. The curve
%  is evaluated at all 31 observable separations; the field is synthesised and
%  rebinned. A gap between them is attributable to the synthesiser rather than
%  to the estimator.


vox = 3; gridsz = [64 64 64];

kmax = 4; V = [];
for dx = 0:kmax, for dy = 0:kmax, for dz = 0:kmax
    if dx == 0 && dy == 0 && dz == 0, continue; end
    V = [V; dx dy dz dx^2+dy^2+dz^2]; %#ok<AGROW>
end, end, end
k2u = unique(V(:,4)); k2u = k2u(k2u <= 48);
REF = readtable(fullfile(OUTV, 'S25_reference_params.csv'));
Kmag = mode_freqs(gridsz, vox);

rows = [];
for w = {'matched','hole'}
    which = w{1};
    C = readtable(fullfile(OUTV, ['S27_ref_' which '_sbcoefs.csv']));
    fprintf('\n===== %s reference =====\n', which);
    fprintf('  %-6s %14s %14s %14s %12s\n', 'shell', 'model gap', 'synthesis gap', ...
            'reference gap', 'range');
    for k = [3 4 5]
        t = REF(strcmp(REF.reference, which) & REF.shell == k, :);
        c = C(C.shell == k, :); p = c.coef; nodes = c.node(2:end);

        g_ref = zeros(numel(k2u),1); g_sb = zeros(numel(k2u),1);
        for j = 1:numel(k2u)
            r = vox*sqrt(k2u(j));
            g_ref(j) = t.nugget(1) + sum(t.sill .* (1 - om3(t.a*r).*exp(-r^2./(2*t.tau.^2))));
            g_sb(j)  = p(1) + sum(p(2:end) .* (1 - om3(nodes*r)));
        end
        S = vg_spectral_field(p, nodes, gridsz, vox);
        g_syn = gen_gamma(S.^2, p(1), V, k2u, gridsz);

        model_gap = max(abs(g_sb  - g_ref));    % estimate's curve vs the truth
        synth_gap = max(abs(g_syn - g_sb));     % what gets simulated vs that curve
        ref_gap   = max(abs(g_syn - g_ref));    % what gets simulated vs the truth
        rng_ = max(g_ref) - min(g_ref);
        fprintf('  %-6d %14.4f %14.4f %14.4f %12.3f\n', k, model_gap, synth_gap, ...
                ref_gap, rng_);
        rows = [rows; strcmp(which,'hole'), k, model_gap, synth_gap, ref_gap, rng_]; %#ok<AGROW>
    end
end

R = array2table(rows, 'VariableNames', ...
    {'is_hole','shell','model_gap','synth_gap','ref_gap','curve_range'});
writetable(R, fullfile(OUTV, 'S30_sb_attribution.csv'));

fprintf('\n================ verdict ================\n');
fprintf('  estimate''s curve against the truth : %.4f at worst (%.2f%% of range)\n', ...
        max(R.model_gap), 100*max(R.model_gap./R.curve_range));
fprintf('  what is simulated against that curve: %.4f at worst (%.2f%% of range)\n', ...
        max(R.synth_gap), 100*max(R.synth_gap./R.curve_range));
if max(R.synth_gap) > 5*max(R.model_gap)
    fprintf('\n  The synthesis loses the agreement, not the model. The error rate in\n');
    fprintf('  S28 is a property of this synthesiser and must not be reported as a\n');
    fprintf('  property of the estimator.\n');
elseif max(R.model_gap) > 5*max(R.synth_gap)
    fprintf('\n  The model loses it. Agreeing on the variogram at 31 separations does\n');
    fprintf('  not imply the same field, and the S28 result stands as a finding.\n');
else
    fprintf('\n  Both contribute at the same order. Neither attribution is safe and\n');
    fprintf('  the S28 admissible-estimate row cannot be interpreted as it stands.\n');
end

%% ========================================================================

% Block 3. Damping the estimate so it can be synthesised
%
%  Replaces each delta shell in the estimate's spectral measure with a damped
%  one,
%
%      C(r) = p0*1{r=0} + sum_j p_j Omega_3(t_j r) exp(-r^2 / (2 tau^2))
%
%  which keeps every term a valid covariance in R^3 and keeps the estimate's
%  own masses and nodes, while spreading each shell over enough lattice modes
%  to be representable. Sweeps tau and reports the departure at each value:
%  large tau preserves the curve, small tau spreads the shells, and the two
%  demands trade off against each other.


vox = 3; gridsz = [64 64 64];

kmax = 4; V = [];
for dx = 0:kmax, for dy = 0:kmax, for dz = 0:kmax
    if dx == 0 && dy == 0 && dz == 0, continue; end
    V = [V; dx dy dz dx^2+dy^2+dz^2]; %#ok<AGROW>
end, end, end
k2u = unique(V(:,4)); k2u = k2u(k2u <= 48);
Kmag = mode_freqs(gridsz, vox);
fprintf('mode spacing %.4f rad/mm, so a bump of width 1/tau needs tau <~ %.0f mm\n', ...
        min(Kmag(Kmag>0)), 1/min(Kmag(Kmag>0)));

REF = readtable(fullfile(OUTV, 'S25_reference_params.csv'));
taus = [10 15 20 30 50 80 200];
rows = [];

for w = {'matched'}
    which = w{1};
    C = readtable(fullfile(OUTV, ['S27_ref_' which '_sbcoefs.csv']));
    for k = [3 4 5]
        t = REF(strcmp(REF.reference, which) & REF.shell == k, :);
        c = C(C.shell == k, :); p = c.coef; nodes = c.node(2:end);
        act = p(2:end) > 1e-12; tk = nodes(act); pk = p(2:end); pk = pk(act);

        g_ref = zeros(numel(k2u),1); g_sb = zeros(numel(k2u),1);
        for j = 1:numel(k2u)
            r = vox*sqrt(k2u(j));
            g_ref(j) = t.nugget(1) + sum(t.sill .* (1 - om3(t.a*r).*exp(-r^2./(2*t.tau.^2))));
            g_sb(j)  = p(1) + sum(pk .* (1 - om3(tk*r)));
        end
        rng_ = max(g_ref) - min(g_ref);
        fprintf('\n===== shell %d, %d active nodes, curve range %.3f =====\n', ...
                k, numel(tk), rng_);
        fprintf('  %8s %14s %16s %14s\n', 'tau', 'curve shift', 'synthesis gap', ...
                'total vs truth');
        for tau = taus
            tv = repmat(tau, numel(tk), 1);
            g_damp = zeros(numel(k2u),1);
            for j = 1:numel(k2u)
                r = vox*sqrt(k2u(j));
                g_damp(j) = p(1) + sum(pk .* (1 - om3(tk*r).*exp(-r^2./(2*tv.^2))));
            end
            S = spec_ref(Kmag, pk, tk, tv, sum(pk));
            g_syn = gen_gamma(S.^2, p(1), V, k2u, gridsz);
            shift = max(abs(g_damp - g_sb));     % what damping costs the curve
            sgap  = max(abs(g_syn  - g_damp));   % what the synthesis still loses
            total = max(abs(g_syn  - g_ref));    % simulated field against the truth
            fprintf('  %8.0f %14.4f %16.4f %14.4f\n', tau, shift, sgap, total);
            rows = [rows; k, tau, shift, sgap, total, rng_]; %#ok<AGROW>
        end
    end
end

R = array2table(rows, 'VariableNames', ...
    {'shell','tau','curve_shift','synth_gap','total_gap','curve_range'});
writetable(R, fullfile(OUTV, 'S31_sb_damped.csv'));

fprintf('\n================ verdict ================\n');
fprintf('  undamped, the simulated field missed the truth by 0.456\n');
best = NaN; bt = NaN;
for tau = taus
    s = R(R.tau == tau, :);
    m = max(s.total_gap);
    fprintf('  tau = %3.0f mm : worst departure from the truth %.4f (%.1f%% of range)\n', ...
            tau, m, 100*max(s.total_gap./s.curve_range));
    if isnan(best) || m < best, best = m; bt = tau; end
end
fprintf('\n  best is tau = %.0f mm at %.4f (%.1f%% of range)\n', bt, best, ...
        100*max(R.total_gap(R.tau == bt)./R.curve_range(R.tau == bt)));
if best < 0.023
    fprintf('  Damping makes the estimate simulable to inside the differences the\n');
    fprintf('  experiment measures. The error-rate claim can be rebuilt on it.\n');
else
    fprintf('  Damping helps but does not reach the accuracy the experiment needs.\n');
    fprintf('  The estimate cannot be simulated faithfully by this route.\n');
end

%% ========================================================================

