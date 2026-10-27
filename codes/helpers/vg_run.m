function [Td, Tc] = vg_run(resdir, fmriprep_root, subject, run, label, ses)
% Depth-shell and fixed-set (cross-eroded) pure-axis variograms for one run.
%
% Returns two tables, both symmetrised over +k and -k:
%   Td  per depth shell   -- tests whether variogram shape varies with distance
%                            from the mask boundary (§6.8)
%   Tc  per erosion radius -- the 6 -> 9 mm decline on a fixed voxel set (§6.7)

SHELLS  = [3 4; 4 5; 5 6; 6 8; 8 99];   % depth from boundary, voxels
KMAX    = 3;                            % lags 3, 6, 9 mm
CROSS_S = [3 4];

if nargin < 5, label = ''; end
if nargin < 6 || isempty(ses), ses = 1; end

[res2, sz, vox, eff] = vg_load(resdir, fmriprep_root, subject, run, ses);
depth = bwdist(~eff);                   % Euclidean distance, voxel units

rows_d = [];
for sh = 1:size(SHELLS,1)
    sel = eff & depth >= SHELLS(sh,1) & depth < SHELLS(sh,2);
    if nnz(sel) < 500, continue; end
    g = vg_axis_gamma(res2, sz, sel, KMAX);
    for a = 1:3
        for k = 1:KMAX
            rows_d = [rows_d; subject, run, SHELLS(sh,1), a, k, ...
                      k*vox(a), g(a,k), nnz(sel)]; %#ok<AGROW>
        end
    end
end

rows_c = [];
for S = CROSS_S
    core = imerode(eff, vg_cross_elem(S));
    if nnz(core) < 500, continue; end
    g = vg_axis_gamma(res2, sz, core, S);
    for a = 1:3
        for k = 1:S
            rows_c = [rows_c; subject, run, S, a, k, ...
                      k*vox(a), g(a,k), nnz(core)]; %#ok<AGROW>
        end
    end
end

Td = array2table(rows_d, 'VariableNames', ...
    {'subject','run','shell_lo','axis','k','distance_mm','gamma','n_vox'});
Tc = array2table(rows_c, 'VariableNames', ...
    {'subject','run','S','axis','k','distance_mm','gamma','n_vox'});
if ~isempty(label)
    Td.model = repmat(string(label), height(Td), 1);
    Tc.model = repmat(string(label), height(Tc), 1);
end
end
