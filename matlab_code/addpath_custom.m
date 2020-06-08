function varargout = addpath_custom (folder, varargin)
%% Add a folder to MATLAB path only if is not already on the path
% Usage: varargout = addpath_custom (folder, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       varargout   - see addpath()
%
% Arguments:
%       folder      - folder to add to path
%                   must be something recognized by addpath()
%       varargin    - see addpath()
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/char2rgb.m
%       cd/plot_grouped_jitter.m
%       cd/plot_violin.m
%       ~/FluoroSNNAP/FluroSNNAP.m

% File History:
% 2019-01-10 Created by Adam Lu
% 2019-11-27 Now checks if the folder exists

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
if iscell(folder)
    [varargout{1:nargout}] = cellfun(@(f) addpath_custom(f, varargin{:}), folder, 'UniformOutput', false);
    return
end

% Add path for each folder if it does not already exist on the MATLAB path
%   and it exists
if ~is_on_path(folder) && isfolder(folder)
    [varargout{1:nargout}] = addpath(folder, varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if ischar(folders)
    folders = {folders};
end

if iscell(folders)
    varargout{:} = ...
        cellfun(@(x) apply_or_return(@(y) is_on_path(y), ...
                                        @(y) addpath(y, varargin{:}), x), ...
                folders, 'UniformOutput', false);
else
    varargout{:} = ...
        arrayfun(@(x) apply_or_return(@(y) is_on_path(y), ...
                                        @(y) addpath(y, varargin{:}), x), ...
                folders);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
