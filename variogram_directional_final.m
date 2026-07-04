%% =====================================================================
%  variogram_directional.m
%  ---------------------------------------------------------------------
%  Computes spatial variograms of fMRI residuals for a set of subjects
%  across multiple runs, then produces both omnidirectional (averaged
%  over runs and sorted by lag distance) and directional outputs.
%
%  PIPELINE
%    1. For each subject and run, load 4D residual images (Res_*.nii)
%       and the corresponding fmriprep brain mask.
%    2. For every lag vector (dx,dy,dz) in one octant up to max_shift,
%       compute the semivariance
%             gamma(h) = 0.5 * mean( (Z(x) - Z(x+h))^2 )
%       using only voxel pairs where both endpoints fall inside the
%       brain mask.
%    3. Average gamma across runs to get the omnidirectional variogram,
%       sort by lag distance in mm, and store per subject.
%    4. For the subjects listed in SUBJECTS_FOR_DIRECTIONAL, also write
%       a long-format CSV preserving (dx, dy, dz, run, gamma) before
%       any averaging or sorting, for downstream anisotropy analysis.
%    5. Re-bin the omnidirectional variograms into nbins uniform
%       distance bins and save those as a separate set of outputs.
%
%  INPUTS
%    rootdir/sub<NN>/run0<R>/Res_*.nii          (SPM residuals)
%    fmriprep_root/sub-<NN>/ses-01/func/        (brain masks)
%
%  OUTPUTS
%    all_subjects_variogram.mat                 raw per-subject gammas
%    all_subjects_variogram.xlsx                same, flattened
%    all_subjects_variogram_binned.mat          uniformly binned
%    all_subjects_variogram_binned.xlsx         same, flattened
%    directional_csvs/sub-<NN>_directional_variogram.csv
%                                               long-format per-shift,
%                                               per-run gammas for the
%                                               subjects in
%                                               SUBJECTS_FOR_DIRECTIONAL
%
%
%  =====================================================================

%% Settings

subjects = [3 5 6 8 9 10 11 17 20 21 26 28 31 32 33 34 36 37 38 39 40];

rootdir       = '<PATH_TO_RESIDUALS_ROOT>';        % contains sub<NN>/run0<R>/Res_*.nii
fmriprep_root = '<PATH_TO_FMRIPREP_DERIVATIVES>';  % contains sub-<NN>/ses-01/func/

runs      = 1:5;
max_shift = 4;   % voxels; max lag ~ max_shift * voxel_size in mm

% Directory for per-subject directional CSVs
DIR_CSV_OUT = fullfile(rootdir, 'directional_csvs');
if ~exist(DIR_CSV_OUT, 'dir'); mkdir(DIR_CSV_OUT); end

% Which subjects to write directional CSVs for. Set to `subjects` to do
% everyone; restricted here to the three representative subjects
% identified by the within-distance dispersion analysis. The
% omnidirectional outputs (.mat / .xlsx) are produced for every subject
% regardless.
SUBJECTS_FOR_DIRECTIONAL = [5 21 28];

all_subjects = struct();

%% Generate shifts (one octant; gamma(h) = gamma(-h))

shifts = [];
for dx = 0:max_shift
    for dy = 0:max_shift
        for dz = 0:max_shift
            if dx == 0 && dy == 0 && dz == 0
                continue;
            end
            shifts = [shifts; dx dy dz];
        end
    end
end

nShifts = size(shifts, 1);

%% Loop over subjects

for si = 1:length(subjects)

    s = subjects(si);
    fprintf('\n============================\n');
    fprintf('Processing Subject %02d\n', s);

    basedir      = fullfile(rootdir, ['sub', sprintf('%02d', s)]);
    fmriprep_dir = fullfile(fmriprep_root, ...
        ['sub-', sprintf('%02d', s)], 'ses-01', 'func');

    all_gammas = zeros(nShifts, length(runs));
    distances  = zeros(nShifts, 1);

    %% Loop over runs

    for r = 1:length(runs)

        run = runs(r);
        fprintf('  Run %d...\n', run);

        resdir = fullfile(basedir, ['run0', num2str(run)]);

        % Load residuals
        files = spm_select('FPList', resdir, '^Res_.*\.nii$');
        V = spm_vol(files);
        res4D = spm_read_vols(V);

        % Load brain mask
        mask_pattern = ['*run-', num2str(run), '*desc-brain_mask.nii.gz'];
        mask_files = dir(fullfile(fmriprep_dir, mask_pattern));

        if isempty(mask_files)
            error(['No mask found for subject ', num2str(s), ', run ', num2str(run)]);
        end

        mask_gz  = fullfile(fmriprep_dir, mask_files(1).name);
        maskfile = mask_gz(1:end-3);

        if ~isfile(maskfile)
            gunzip(mask_gz);
        end

        mask = niftiread(maskfile) > 0;

        % Check dimensions
        if ~isequal(size(mask), size(res4D(:,:,:,1)))
            error('Mask and residual dimensions do not match');
        end

        % Voxel size (mm) from affine
        voxel_size = sqrt(sum(V(1).mat(1:3,1:3).^2));
        T = size(res4D, 4);

        % Variogram
        for i = 1:nShifts

            h = shifts(i, :);

            shifted      = circshift(res4D, [h 0]);
            shifted_mask = circshift(mask, h);

            valid   = mask & shifted_mask;
            valid4D = repmat(valid, [1 1 1 T]);

            diff_sq = (res4D - shifted).^2;
            diff_sq = diff_sq(valid4D);
            diff_sq = diff_sq(~isnan(diff_sq));
            n_pairs(i, r) = numel(diff_sq);

            all_gammas(i, r) = 0.5 * mean(diff_sq);

            if r == 1
                distances(i) = sqrt(sum((h .* voxel_size).^2));
            end
        end
    end

    %% Average across runs

    mean_gamma = mean(all_gammas, 2);
    std_gamma  = std(all_gammas, 0, 2);

    % Sort by lag distance
    [distances_sorted, idx] = sort(distances);
    mean_gamma = mean_gamma(idx);
    std_gamma  = std_gamma(idx);

    % Store
    all_subjects(si).subject    = s;
    all_subjects(si).distances  = distances_sorted;
    all_subjects(si).mean_gamma = mean_gamma;
    all_subjects(si).std_gamma  = std_gamma;

    %% Write directional CSV for selected subjects

    if ismember(s, SUBJECTS_FOR_DIRECTIONAL)

        % Long-format table: one row per (shift, run).
        nRows = nShifts * length(runs);
        Subject     = repmat(s, nRows, 1);
        run_col     = zeros(nRows, 1);
        dx_col      = zeros(nRows, 1);
        dy_col      = zeros(nRows, 1);
        dz_col      = zeros(nRows, 1);
        distance_mm = zeros(nRows, 1);
        gamma_run   = zeros(nRows, 1);

        row = 1;
        for r = 1:length(runs)
            for i = 1:nShifts
                run_col(row)     = runs(r);
                dx_col(row)      = shifts(i, 1);
                dy_col(row)      = shifts(i, 2);
                dz_col(row)      = shifts(i, 3);
                distance_mm(row) = sqrt(sum((shifts(i,:) .* voxel_size).^2));
                gamma_run(row)   = all_gammas(i, r);
                row = row + 1;
            end
        end

        dir_table = table(Subject, run_col, dx_col, dy_col, dz_col, ...
                          distance_mm, gamma_run, ...
                          'VariableNames', ...
                          {'Subject','run','dx','dy','dz','distance_mm','gamma_run'});

        out_csv = fullfile(DIR_CSV_OUT, ...
                           sprintf('sub-%02d_directional_variogram.csv', s));
        writetable(dir_table, out_csv);
        fprintf('  Wrote directional CSV: %s\n', out_csv);
    end

end

%% Grid plot of raw variograms

nSub  = length(subjects);
ncols = 5;
nrows = ceil(nSub / ncols);

figure;
for i = 1:nSub
    subplot(nrows, ncols, i);
    errorbar(all_subjects(i).distances, ...
             all_subjects(i).mean_gamma, ...
             all_subjects(i).std_gamma, ...
             'o-', 'LineWidth', 1);
    title(['Sub ', num2str(all_subjects(i).subject)]);
    xlabel('Lag distance (mm)');
    ylabel('Semivariance \gamma(h)');
    grid on;
end
sgtitle('Variograms for All Subjects');

%% Save raw outputs

save('all_subjects_variogram.mat', 'all_subjects');

excel_data = table();
for i = 1:nSub
    s  = all_subjects(i).subject;
    d  = all_subjects(i).distances;
    g  = all_subjects(i).mean_gamma;
    sd = all_subjects(i).std_gamma;
    n  = length(d);

    omni_table = table( ...
        repmat(s, n, 1), d, g, sd, ...
        'VariableNames', {'Subject','Distance','MeanGamma','StdGamma'});
    excel_data = [excel_data; omni_table];
end
writetable(excel_data, 'all_subjects_variogram.xlsx');

fprintf('\nDone: all subjects processed. Directional CSVs in %s\n', DIR_CSV_OUT);

%% Uniform binning

load('all_subjects_variogram.mat');  % loads "all_subjects"

nbins = 20;   % uniform-width distance bins for the binned variogram

binned_subjects = struct();

for i = 1:length(all_subjects)

    s         = all_subjects(i).subject;
    distances = all_subjects(i).distances(:);
    gamma     = all_subjects(i).mean_gamma(:);

    % Create uniform bins
    edges       = linspace(min(distances), max(distances), nbins + 1);
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

    binned_gamma = nan(nbins, 1);
    binned_std   = nan(nbins, 1);
    counts       = zeros(nbins, 1);

    % Binning loop
    for b = 1:nbins
        idx = distances >= edges(b) & distances < edges(b+1);
        if any(idx)
            vals             = gamma(idx);
            binned_gamma(b)  = mean(vals);
            binned_std(b)    = std(vals);
            counts(b)        = sum(idx);
        end
    end

    % Remove empty bins
    valid        = ~isnan(binned_gamma);
    bin_centers  = bin_centers(valid);
    binned_gamma = binned_gamma(valid);
    binned_std   = binned_std(valid);
    counts       = counts(valid);

    % Safety alignment (force column vectors of equal length)
    bin_centers  = bin_centers(:);
    binned_gamma = binned_gamma(:);
    binned_std   = binned_std(:);
    counts       = counts(:);

    n_bins_valid = min([length(bin_centers), length(binned_gamma), ...
                        length(binned_std), length(counts)]);

    bin_centers  = bin_centers(1:n_bins_valid);
    binned_gamma = binned_gamma(1:n_bins_valid);
    binned_std   = binned_std(1:n_bins_valid);
    counts       = counts(1:n_bins_valid);

    % Store
    binned_subjects(i).subject     = s;
    binned_subjects(i).bin_centers = bin_centers;
    binned_subjects(i).mean_gamma  = binned_gamma;
    binned_subjects(i).std_gamma   = binned_std;
    binned_subjects(i).counts      = counts;

end

%% Grid plot of binned variograms

nSub  = length(binned_subjects);
ncols = 5;
nrows = ceil(nSub / ncols);

figure;
for i = 1:nSub
    subplot(nrows, ncols, i);
    errorbar(binned_subjects(i).bin_centers, ...
             binned_subjects(i).mean_gamma, ...
             binned_subjects(i).std_gamma, ...
             'o-', 'LineWidth', 1);
    title(['Sub ', num2str(binned_subjects(i).subject)]);
    xlabel('Lag distance (mm)');
    ylabel('Semivariance \gamma(h)');
    grid on;
end
sgtitle('Binned Variograms (Uniform Binning)');

%% Save binned outputs

save('all_subjects_variogram_binned.mat', 'binned_subjects');

excel_data = table();
for i = 1:nSub
    s  = binned_subjects(i).subject;
    d  = binned_subjects(i).bin_centers(:);
    g  = binned_subjects(i).mean_gamma(:);
    sd = binned_subjects(i).std_gamma(:);
    c  = binned_subjects(i).counts(:);

    % Ensure equal length (final safeguard)
    n_rows = min([length(d), length(g), length(sd), length(c)]);
    d  = d(1:n_rows);
    g  = g(1:n_rows);
    sd = sd(1:n_rows);
    c  = c(1:n_rows);

    binned_table = table( ...
        repmat(s, n_rows, 1), d, g, sd, c, ...
        'VariableNames', {'Subject','Distance','MeanGamma','StdGamma','Counts'});
    excel_data = [excel_data; binned_table];
end
writetable(excel_data, 'all_subjects_variogram_binned.xlsx');

fprintf('\nDone: binning complete, files saved successfully.\n');
