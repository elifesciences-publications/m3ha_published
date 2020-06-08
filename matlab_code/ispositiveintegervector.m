function isPositiveIntegerVector = ispositiveintegervector (x)
%% Returns whether an input is a positive integer vector
% Usage: isPositiveIntegerVector = ispositiveintegervector (x)
% Explanation:
%       Tests whether the input is a positive integer vector
% Example(s):
%       ispositiveintegervector(1:10)
%       ispositiveintegervector(magic(3))
%       ispositiveintegervector(-1:3)
% Outputs:
%       isPositiveIntegerVector
%                       - whether the input is a positive integer vector
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires: 
%       cd/create_error_for_nargin.m
%       cd/isaninteger.m
%
% Used by:
%       cd/extract_param_values.m
%       cd/extract_vars.m
%       cd/match_dimensions.m
%       cd/match_time_info.m
%       cd/m3ha_select_cells.m
%       cd/m3ha_plot_bar3.m
%       cd/m3ha_plot_violin.m
%       cd/plot_struct.m
%       cd/run_neuron.m
%       cd/set_axes_properties.m

% File History:
% 2018-10-24 Created by Adam Lu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isPositiveIntegerVector = isnumeric(x) && isvector(x) && ...
                            all(isaninteger(x) & x > 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%