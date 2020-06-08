function array = vec2array (vector, dims)
%% Convert a vector to an array with dimensions given by dims using linear indexing (obsolete, use reshape() instead)
% Usage: array = vec2array (vector, dims)
% Outputs:  
%       array       - array with dimensions given by dims
% Arguments:    
%       vector      - a vector
%                   must be a vector
%       dims        - dimensions for the output array
%                   must be a vector of positive integers
%
% Used by:
%       cd/m3ha_network_tuning_maps.m
%       /media/adamX/RTCl/tuning_maps.m

% 2017-05-04 Created by Adam Lu
% 2017-05-05 Expanded to accept cell arrays and struct arrays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an input Parser
addRequired(iP, 'vector', @isvector);       % a vector
addRequired(iP, 'dims', ...                 % dimensions for the output array
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer'}));

% Read from the input Parser
parse(iP, vector, dims);

% Check relationships between arguments
if (prod(dims) ~= length(vector))
    error('product of dimensions must equal length of vector!');
end

%% Initialize array
if isnumeric(vector)
    array = zeros(dims);
elseif iscell(vector)
    array = cell(dims);
elseif isstruct(vector)
    array = repmat(struct, dims);
end

%% Map each element of the vector to an element in the array 
%   with built-in Matlab linear indexing
for k = 1:length(vector)
    if isnumeric(vector)
        array(ind2sub(dims, k)) = vector(k);
    elseif iscell(vector) || isstruct(vector)
        array{ind2sub(dims, k)} = vector{k};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%