%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate combinations with repetation
% Solution provided by Eitan T: http://stackoverflow.com/questions/18591440/how-to-find-all-permutations-with-repetition-in-matlab
% simplyfing Amro's anwser: http://stackoverflow.com/questions/4165859/matlab-generate-all-possible-combinations-of-the-elements-of-some-vectors
% Input: v is set of possible letters and k is Length of each permutation
% Output: All possible  combinations with repetation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[X]=combrep(v,k)

% Create all possible permutations (with repetition) of letters stored in v
C = cell(k, 1);             %// Preallocate a cell array
[C{:}] = ndgrid(v);         %// Create K grids of values
X = cellfun(@(x){x(:)}, C); %// Convert grids to column vectors
X = [X{:}];                 %// Obtain all permutations

end