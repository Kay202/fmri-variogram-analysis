function [res2, sz, vox, eff] = vg_load(resdir, fmriprep_root, subject, run, ses)
% Load SPM residuals and the effective analysis mask for one run.
%
%   res2  (T x nvox) single, NaNs zeroed
%   sz    [X Y Z]
%   vox   voxel size in mm
%   eff   logical volume: fMRIPrep brain mask AND finite residuals
%   ses   optional session number for the mask lookup (default 1)
%
% Accepts either residual layout:
%   * SPM's native 482 float64 3D volumes, Res_*.nii   (the F4 analyses)
%   * one float32 4D file, Res4D_*.nii                 (residuals_modelB)
% The two are numerically identical at float32 precision -- verified to a max
% relative difference of 5.9e-8, below float32 eps.
%
% Time is kept down the rows because MATLAB is column-major: res2(:,idx) then
% gathers contiguous memory, which is what makes the per-lag loop tolerable.

if nargin < 5 || isempty(ses), ses = 1; end

f4d = dir(fullfile(resdir, 'Res4D_*.nii'));
if ~isempty(f4d)
    fname = fullfile(resdir, f4d(1).name);
    V = spm_vol(fname);                 % one struct per volume of the 4D file
    Y = spm_read_vols(V);
else
    files = spm_select('FPList', resdir, '^Res_.*\.nii$');
    assert(~isempty(files), 'no Res_*.nii or Res4D_*.nii in %s', resdir);
    V = spm_vol(files);
    Y = spm_read_vols(V);
end

d  = size(Y);
T  = d(4);
sz = d(1:3);
vox = sqrt(sum(V(1).mat(1:3,1:3).^2));

res2 = single(reshape(Y, [], T))';
clear Y

finite_mask = reshape(all(isfinite(res2), 1), sz);
res2(~isfinite(res2)) = 0;

fd = fullfile(fmriprep_root, sprintf('sub-%02d', subject), sprintf('ses-%02d', ses), 'func');
f  = dir(fullfile(fd, sprintf('*run-%d*desc-brain_mask.nii', run)));
if isempty(f)
    gz = dir(fullfile(fd, sprintf('*run-%d*desc-brain_mask.nii.gz', run)));
    assert(~isempty(gz), 'no brain mask for sub-%02d run %d', subject, run);
    gunzip(fullfile(fd, gz(1).name));
    f = dir(fullfile(fd, sprintf('*run-%d*desc-brain_mask.nii', run)));
end
bmask = niftiread(fullfile(fd, f(1).name)) > 0;
assert(isequal(size(bmask), sz), 'mask/residual dims differ');

eff = bmask & finite_mask;
end
