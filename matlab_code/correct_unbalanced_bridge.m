function vvecNew = correct_unbalanced_bridge (vvecOld, ivecOld, varargin)
%% Shifts a current pulse response to correct the unbalanced bridge
% Usage: vvecNew = correct_unbalanced_bridge (vvecOld, ivecOld, varargin)
% Arguments:    
%       vvecOld     - Voltage vector to correct
%                   must be a numeric vector
%       ivecOld     - Current vector
%                   must be a numeric vector
%       varargin    - 'UseCurrentFlag': whether to use the current trace
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/compute_initial_slopes.m
%
% Used by:
%       cd/m3ha_correct_unbalanced_bridge.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_initial_slopes.m

% File History:
% 2018-07-25 BT - Created by Brian
% 2018-08-10 AL - Now checks number of arguments
% 2018-08-11 AL - Added UseCurrentFlag and set the default to not use it
% 2018-08-12 AL - Set default of UseCurrentFlag to true
% 2018-08-13 AL - Now uses the actual endpoints used by 
%                   compute_initial_slopes to find the indices to shift

%% Hard-coded parameters
nSamples = 2;       % number of samples between the voltage jump
                    %   caused by an unbalanced bridge

%% Default values for optional arguments
useCurrentFlagDefault = true;      % use the current trace by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vvecOld', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'ivecOld', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'UseCurrentFlag', useCurrentFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vvecOld, ivecOld, varargin{:});
useCurrentFlag = iP.Results.UseCurrentFlag;

%% Do the job
% Initialize the new voltage vector as the old one
vvecNew = vvecOld;

% Create an artificial time vector
tvecOld = 1:length(ivecOld);

% Get the average initial slope
if useCurrentFlag
    [avgSlope, ~, ~, indsUsed] = ...
        compute_initial_slopes(tvecOld, vvecOld, ...
                               'IvecCpr', ivecOld, ...
                               'NSamples', nSamples);
else
    [avgSlope, ~, ~, indsUsed] = ...
        compute_initial_slopes(tvecOld, vvecOld, ...
                               'NSamples', nSamples);
end

% Return if there is no pulse or average slope computed
if isnan(avgSlope)
    fprintf('There is no pulse or average slope computed!!\n');
    return
end

% Get the indices for the region to shift
idxLast1 = indsUsed(2);
idxFirst2 = indsUsed(3);
indToShift = idxLast1:idxFirst2;

% Calculate the time difference that the average slope was computed for
deltaTime = tvecOld(nSamples) - tvecOld(1);

% Calculate the voltage to correct by
shiftBy = avgSlope * deltaTime;

% Correct voltage
vvecNew(indToShift) = vvecOld(indToShift) + shiftBy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

function [vvecNew] = correct_unbalanced_bridge (tvecOld, vvecOld, ivecOld, nSamples)

parse(iP, tvecOld, vvecOld, ivecOld, nSamples);

addRequired(iP, 'tvecOld', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'nSamples', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Get index of current pulse endpoints
[idxStart, idxEnd] = find_pulse_endpoints(ivecOld);

shiftBy = avgSlope * (tvecOld(idxStart + nSamples - 1) - tvecOld(idxStart));

vvecNew(idxStart + 1:idxEnd) = vvecOld(idxStart + 1:idxEnd) + shiftBy;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%