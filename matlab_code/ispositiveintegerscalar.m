function isPositiveIntegerScalar = ispositiveintegerscalar (x)
%% Returns whether an input is a positive integer scalar
% Usage: isPositiveIntegerScalar = ispositiveintegerscalar (x)
% Explanation:
%       Tests whether the input is a positive integer scalar
% Example(s):
%       ispositiveintegerscalar(10)
%       ispositiveintegerscalar(10.5)
%       ispositiveintegerscalar(-1:3)
% Outputs:
%       isPositiveIntegerScalar
%                       - whether the input is a positive integer scalar
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires: 
%       cd/isaninteger.m
%
% Used by:
%       cd/compute_weighted_average.m
%       cd/convert_to_char.m
%       cd/detect_spikes_multiunit.m
%       cd/find_custom.m
%       cd/find_in_strings.m
%       cd/m3ha_neuron_create_sim_params.m
%       cd/match_positions.m
%       cd/plot_fitted_traces.m
%       cd/plot_grouped_histogram.m
%       cd/plot_grouped_scatter.m
%       cd/plot_struct.m
%       cd/plot_traces.m
%       cd/read_line_from_file.m
%       cd/renamevars_custom.m
%       ~/m3ha/optimizer4gabab/singleneuronfitting42.m and beyond

% File History:
% 2018-10-26 Adapted from ispositiveintegervector.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isPositiveIntegerScalar = isnumeric(x) && isscalar(x) && ...
                            isaninteger(x) && x > 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
