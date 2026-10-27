function p = vg_paths(varargin)
%VG_PATHS  Every filesystem location the imaging stage needs.
%
%   p = vg_paths()            all paths, validated
%   p = vg_paths('spm')       one path, validated
%   p = vg_paths('spm','logs')
%
% Three locations depend on your machine and cannot be guessed. Set them as
% environment variables, or edit the defaults below.
%
%   VG_SPM        SPM installation (SPM 25.01.02 was used)
%   VG_FMRIPREP   fMRIPrep derivatives, the directory holding sub-*/
%   VG_LOGS       Presentation logs, the directory holding Sub*/ExpInfo_run*.mat
%   VG_RESIDUALS  where to put the Model B residual store. About 132 GB, so
%                 it belongs outside this package.
%
% Two more are optional and only needed by scripts that read the original
% first-level model rather than the Model B residuals:
%
%   VG_SPMFIRST   SPM first-level output, the directory holding sub*/run*/
%
% Everything else is derived from the project root, which is taken to be two
% directories above this file and can be overridden with VG_ROOT:
%
%   p.root        project root
%   p.residuals   Model B residual store, written by GLM.m
%   p.out         out_matlab, where the variogram stage writes its CSVs
%
% A path that is needed and not set raises an error naming the variable, rather
% than failing later on a missing file.
% edit here if you would rather not set environment variables
defaults.spm      = '';
defaults.fmriprep = '';
defaults.logs     = '';
defaults.spmfirst = '';
defaults.residuals = '';
here = fileparts(mfilename('fullpath'));
root = getenv('VG_ROOT');
if isempty(root), root = fileparts(fileparts(here)); end
p.root      = root;
p.residuals = pick('VG_RESIDUALS', defaults.residuals);
p.out       = fullfile(root, 'data', 'out_matlab');
p.helpers   = here;
p.spm      = pick('VG_SPM',      defaults.spm);
p.fmriprep = pick('VG_FMRIPREP', defaults.fmriprep);
p.logs     = pick('VG_LOGS',     defaults.logs);
p.spmfirst = pick('VG_SPMFIRST', defaults.spmfirst);
needed = varargin;
if isempty(needed), needed = {'spm', 'fmriprep', 'logs'}; end
describe = struct( ...
    'spm',      'SPM installation (VG_SPM)', ...
    'fmriprep', 'fMRIPrep derivatives (VG_FMRIPREP)', ...
    'logs',     'Presentation logs (VG_LOGS)', ...
    'spmfirst', 'SPM first-level output (VG_SPMFIRST)', ...
    'residuals','Model B residual store (VG_RESIDUALS)', ...
    'out',      'out_matlab output directory', ...
    'root',     'project root');
for i = 1:numel(needed)
    k = needed{i};
    if ~isfield(p, k)
        error('vg_paths: no such path "%s"', k);
    end
    if isempty(p.(k))
        error(['vg_paths: %s is not set.\n' ...
               'Set the environment variable or edit the defaults at the top ' ...
               'of %s.'], describe.(k), [mfilename '.m']);
    end
    if ~strcmp(k, 'out') && exist(p.(k), 'dir') ~= 7
        error('vg_paths: %s is set to "%s", which is not a directory.', ...
              describe.(k), p.(k));
    end
end
if ~isempty(p.out) && exist(p.out, 'dir') ~= 7, mkdir(p.out); end
end
function v = pick(envvar, dflt)
v = getenv(envvar);
if isempty(v), v = dflt; end
end
