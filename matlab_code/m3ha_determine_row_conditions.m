function rowConditions = ...
                m3ha_determine_row_conditions (rowMode, pharmToUse, gIncrToUse)
%% Determine the conditions for each row
% Usage: rowConditions = ...
%               m3ha_determine_row_conditions (rowMode, pharmToUse, gIncrToUse)
% 
% Explanation: 
%       TODO
% Example(s):
%       TODO
% Outputs:
%       TODO
%
% Arguments:
%       TODO
%
% Requires:
%
% Used by:
%       cd/m3ha_select_raw_traces.m

% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-05-22 Changed line width and indentation
% 2018-11-15 Moved to Adams_Functions
% 2018-11-15 Improved documentation
% 2019-12-21 No longer depends on colMode and attemptNumber
% TODO: Make this depend only on nPCond and nGIncr?

%% Hard-coded parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Count the number of possible pharm conditions 
nPCond = length(pharmToUse);

% Count the number of possible conductance amplitude scaling percentages
nGIncr = length(gIncrToUse);

%% Do the job
if rowMode == 1
    % Each row is one pharm condition
    rowConditions = pharmToUse;
elseif rowMode == 2
    % Each row is a pharm condition paired with a g incr
    nRows = nPCond * nGIncr;

    % Find each row condition
    rowConditions = zeros(nRows, 2);
    for iRow = 1:nRows
        % Find pharm condition of each row (index in pharmToUse)
        rowConditions(iRow, 1) = floor((iRow-1) / nGIncr) + 1;    

        % Find g-incr condition of each row (index in gIncrToUse)
        rowConditions(iRow, 2) = mod(iRow - 1, nGIncr) + 1;        
    end
else
    error('row mode undefined!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%