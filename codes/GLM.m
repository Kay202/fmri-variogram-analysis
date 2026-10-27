function GLM(subjects, runs)
% Model B first level for both sessions, saving residuals locally.
%
% Model B is the extended-nuisance specification:
%   6 task regressors (5 deviants + adaptor train)
%   6 motion + 6 motion derivatives + 6 aCompCor        = 18 nuisance
%   AR(1), hpf 128, mthresh 0.8, no explicit mask, global None
%
% Run index r = 1..10 follows the authors' convention: r <= 5 is ses-01 run r,
% r > 5 is ses-02 run r-5. Conditions are rebuilt from the original
% Presentation logs for every run, since the release has no events.tsv files.
%
% SPM writes residuals as 482 separate float64 volumes per run, 1.25 GB, or
% 263 GB over the 210 runs. Each run is therefore converted to a single
% float32 4D file of 627 MB, about 132 GB in total. Only residuals are kept;
% betas, ResMS and RPV are deleted.

if nargin < 1 || isempty(subjects)
    subjects = [3 5 6 8 9 10 11 17 20 21 26 28 31 32 33 34 36 37 38 39 40];
end
if nargin < 2 || isempty(runs), runs = 1:10; end

addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
P        = vg_paths('spm', 'logs', 'fmriprep');
SPMDIR   = P.spm;
LOGROOT  = P.logs;
FMRIPREP = P.fmriprep;
OUT      = P.residuals;
MIN_FREE_GB = 20;   % refuse to continue if the disk gets this tight

addpath(SPMDIR);
if ~exist(OUT,'dir'), mkdir(OUT); end
spm('defaults','FMRI'); spm_jobman('initcfg');

condcode  = (1:7)';
condnames = {'identical';'superclose';'close';'distant';'superdistant'; ...
             'fixation';'adaptation'};
condtable = table(condcode, condnames);
interest  = 1:5;

t0all = tic; nDone = 0; failed = {};
for s = subjects
    for r = runs
      try
        if r <= 5, ss = 1; rr = r; else, ss = 2; rr = r - 5; end
        od = fullfile(OUT, sprintf('sub%02d', s), sprintf('run%02d', r));
        out4d = fullfile(od, sprintf('Res4D_sub%02d_run%02d.nii', s, r));

        if isfile(out4d)
            fprintf('sub-%02d run%02d: already done, skipping\n', s, r); continue;
        end

        free = disk_free_gb(OUT);
        if free < MIN_FREE_GB
            error('only %.1f GB free -- stopping before the disk fills', free);
        end

        fprintf('\n===== sub-%02d ses-%02d run-%d (log run%02d) | %.0f GB free =====\n', ...
                s, ss, rr, r, free);
        t0 = tic;
        if ~exist(od,'dir'), mkdir(od); end

        % conditions from the original Presentation log
        L = load(fullfile(LOGROOT, ['Sub' num2str(s)], sprintf('ExpInfo_run%02d.mat', r)));
        [names, onsets, durations] = build_conditions(L.ExpInfo, condtable, interest); %#ok<ASGLU>
        condpath = fullfile(od,'conditions.mat');
        save(condpath, 'names', 'onsets', 'durations');

        % model B nuisance
        conf = fullfile(FMRIPREP, sprintf('sub-%02d',s), sprintf('ses-%02d',ss), 'func', ...
            sprintf(['sub-%02d_ses-%02d_task-releasesuppression_run-%d_' ...
                     'desc-confounds_timeseries.tsv'], s, ss, rr));
        tb = readtable(conf,'FileType','text','Delimiter','\t', ...
                       'TreatAsMissing',{'n/a','NA',''});
        R = [tb.trans_x tb.trans_y tb.trans_z tb.rot_x tb.rot_y tb.rot_z, ...
             tb.trans_x_derivative1 tb.trans_y_derivative1 tb.trans_z_derivative1, ...
             tb.rot_x_derivative1   tb.rot_y_derivative1   tb.rot_z_derivative1, ...
             tb.a_comp_cor_00 tb.a_comp_cor_01 tb.a_comp_cor_02, ...
             tb.a_comp_cor_03 tb.a_comp_cor_04 tb.a_comp_cor_05];
        R(~isfinite(R)) = 0;
        regpath = fullfile(od,'regressors.txt');
        writematrix(R, regpath, 'Delimiter','tab');

        % specify + estimate
        bold = fullfile(FMRIPREP, sprintf('sub-%02d',s), sprintf('ses-%02d',ss), 'func', ...
            sprintf(['sub-%02d_ses-%02d_task-releasesuppression_run-%d_' ...
                     'space-MNI152NLin2009cAsym_desc-preproc_bold.nii'], s, ss, rr));
        if ~isfile(bold)
            gz = [bold '.gz'];
            if ~isfile(gz)
                fprintf('  MISSING bold, skipped\n');
                failed{end+1} = sprintf('sub-%02d run%02d (no bold)', s, r); %#ok<AGROW>
                continue;
            end
            % gunzip reads the placeholder directly, so hydrate it first --
            % an unguarded gunzip on a cold .gz killed an earlier run outright
            if ~hydrate_file(gz)
                fprintf('  could not hydrate .gz, skipped\n');
                failed{end+1} = sprintf('sub-%02d run%02d (gz hydration)', s, r); %#ok<AGROW>
                continue;
            end
            gunzip(gz);
        end
        scans = cellstr(hydrate_and_expand(bold));
        if isfile(fullfile(od,'SPM.mat')), delete(fullfile(od,'SPM.mat')); end

        clear matlabbatch
        matlabbatch{1}.spm.stats.fmri_spec.dir            = {od};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units   = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = 1;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans     = scans;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond      = ...
            struct('name',{},'onset',{},'duration',{},'tmod',{},'pmod',{},'orth',{});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi     = {condpath};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress   = struct('name',{},'val',{});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {regpath};
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf       = 128;
        matlabbatch{1}.spm.stats.fmri_spec.fact           = struct('name',{},'levels',{});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt           = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global         = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh        = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask           = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi            = 'AR(1)';
        matlabbatch{1}.spm.stats.fmri_spec.write_residuals = 1;
        spm_jobman('run', matlabbatch);

        clear matlabbatch
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(od,'SPM.mat')};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('run', matlabbatch);

        % 482 float64 3D volumes -> one float32 4D file
        nv = pack_residuals(od, out4d);

        % keep residuals only
        delete(fullfile(od,'Res_*.nii'));
        delete(fullfile(od,'beta_*.nii'));
        if isfile(fullfile(od,'ResMS.nii')), delete(fullfile(od,'ResMS.nii')); end
        if isfile(fullfile(od,'RPV.nii')),   delete(fullfile(od,'RPV.nii'));   end

        nDone = nDone + 1;
        fprintf('  packed %d volumes -> %.0f MB | %.1f min | %d runs done, %.1f h elapsed\n', ...
                nv, dir(out4d).bytes/1e6, toc(t0)/60, nDone, toc(t0all)/3600);

      catch ME
        % one transient OneDrive read should not end a six-hour job; log the
        % run, leave it unpacked so a rerun picks it up, and carry on
        if strcmp(ME.identifier, 'MATLAB:nomem') || contains(ME.message, 'GB free')
            rethrow(ME);          % genuine resource limits still stop us
        end
        fprintf('  !! sub-%02d run%02d FAILED: %s\n', s, r, ME.message);
        failed{end+1} = sprintf('sub-%02d run%02d (%s)', s, r, ME.identifier); %#ok<AGROW>
      end
    end
end
fprintf('\nDONE: %d runs in %.1f h\n', nDone, toc(t0all)/3600);
if isempty(failed)
    fprintf('no failures\n');
else
    fprintf('%d runs need a rerun:\n', numel(failed));
    fprintf('  %s\n', failed{:});
end
end

function [names, onsets, durations] = build_conditions(ExpInfo, condtable, interest)
fr = ExpInfo.Cfg.Screen.FrameRate;
n = numel(ExpInfo.TrialInfo);
trialcode = zeros(n,1); trialstart = zeros(n,1);
adaptdur = cell(n,1); deviantdur = cell(n,1);
for t = 1:n
    trialcode(t)  = ExpInfo.TrialInfo(t).trial.code;
    trialstart(t) = ExpInfo.TrialInfo(t).tStart;
    if trialcode(t) ~= 6
        adaptdur{t}   = ExpInfo.TrialInfo(t).trial.pageDuration(1:end-4);
        deviantdur{t} = ExpInfo.TrialInfo(t).trial.pageDuration(end-3:end);
    end
end
names = cell(0); onsets = cell(0); durations = cell(0);
aon = []; adu = [];
for code = interest
    tr = find(trialcode == code);
    cn = find(condtable.condcode == code);
    a0 = zeros(numel(tr),1); a1 = a0; d0 = a0; d1 = a0;
    for m = 1:numel(tr)
        a0(m) = trialstart(tr(m));
        a1(m) = sum(adaptdur{tr(m)})/fr;
        d0(m) = a0(m) + a1(m);
        d1(m) = sum(deviantdur{tr(m)})/fr;
    end
    aon = [aon; a0]; adu = [adu; a1]; %#ok<AGROW>
    names = [names; condtable.condnames(cn)]; %#ok<AGROW>
    onsets = [onsets; {d0}]; durations = [durations; {d1}]; %#ok<AGROW>
end
[aon, ix] = sort(aon); adu = adu(ix);
names = [names; condtable.condnames(7)];
onsets = [onsets; {aon}]; durations = [durations; {adu}];
end

function nv = pack_residuals(od, out4d)
% SPM writes float64 3D volumes; store one float32 4D file instead. Halves the
% bytes and takes 482 files to 1, at no cost to the variogram (which computes in
% float32 and reproduces the published curve to 2e-8).
%
% Written with a single bulk niftiwrite rather than 482 spm_write_vol calls --
% the per-volume loop took ~8.5 min per run and dominated the whole job.
files = spm_select('FPList', od, '^Res_.*\.nii$');
V = spm_vol(files);
nv = numel(V);
Y = single(spm_read_vols(V));

info = niftiinfo(deblank(files(1,:)));
info.ImageSize       = size(Y);
info.PixelDimensions = [info.PixelDimensions(1:3) 1];
info.Datatype        = 'single';
info.BitsPerPixel    = 32;
info.raw.dim         = [4 size(Y) 1 1 1];
info.raw.datatype    = 16;
info.raw.bitpix      = 32;
info.Description     = 'model B residuals, float32 4D';

[p, n, e] = fileparts(out4d);
if strcmpi(e, '.nii'), stem = fullfile(p, n); else, stem = out4d; end
niftiwrite(Y, stem, info, 'Compressed', false);
end

function g = disk_free_gb(p)
d = java.io.File(p);
g = d.getFreeSpace() / 1e9;
end

function ok = hydrate_file(f)
% Force a OneDrive placeholder local by touching its first bytes, retrying a few
% times. Cold reads on this share fail outright often enough to matter.
ok = false;
for attempt = 1:5
    try
        fid = fopen(f, 'r');
        if fid > 0
            n = numel(fread(fid, 4096, 'uint8'));
            fclose(fid);
            if n > 0, ok = true; return; end
        end
    catch
    end
    fprintf('  hydration attempt %d on %s failed, retrying...\n', attempt, ...
            char(java.io.File(f).getName()));
    pause(20);
end
end

function files = hydrate_and_expand(bold)
for attempt = 1:5
    try
        fid = fopen(bold, 'r');
        if fid > 0, fread(fid, 4096, 'uint8'); fclose(fid); end
        files = spm_select('expand', bold);
        if size(files,1) > 1, return; end
    catch ME
        fprintf('  hydration attempt %d failed (%s), retrying...\n', attempt, ME.identifier);
    end
    pause(20);
end
error('could not read %s after 5 attempts', bold);
end
