function  varargout = handleVarargin(optionalInputs, optionalVariables, defaultValues)
% assignedInputs = assignDefaultValues(optionalInputs, optionalVariables, defaultValues)
% This function can be used to when you 'varargin'. It will assing
% varargin{1, 2, ...} to their corresponding variables (you don't need to
% do if nargin > 0 etc). If some of the optional variables are not passed
% then the default values will be assigned.
% All the optinal variables (optionalVariables) should be bundlled in a
% structure (which are assigend to empty before be passed to the function).
% All the default values (defaultValues) are bundled in a cell array. 
% All the 'varargin' already bundelled in a cell which shlould be passed to
% the function as the 'optionalInputs'. 
% Those optional inputs (optionalInputs) which are not passed or passed
% empty will be assigned by default values.
% 
% EXAMPLE:
% Considering dealing with this function:
% generateHayashiStim(imageName, phaseShiftValue, varargin)
% in which there are 3 optinal input: imagePath, outputPath and movieName
% generateHayashiStim(imageName, phaseShiftValue, imagePath, outputPath, movieName)
% 
% To assign the optinal/default values we need to provide the variable
% names:
% optionalVariables.imagePath = []; optionalVariables.outputPath = []; optionalVariables.movieName = []; 
% and their default values:
% defaultValues{1} = '.'; defaultValues{2} = 'B\folder1\subfolder1'; defaultValues{3} = 'movieName';
% and pass the to the function
% optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);
% 
% ------
% Input:
% 1     optionalInputs: cell array
%       It's your varargin cell array
% 2     optionalVariables: structure
%       all the optinal variables (preassigined to empty) bundled in a
%       structure
% 3     defaultValues: cell array
%       default values with the same order with 'optional variables' and 'optional inputs'
%       
% Output: 
% (1)	assignedInputs: structure
%       optional varibales (bundles in output structure) with default
%       values according to user/default values
%
% ------
% potential improvments:
% (1) check if the number of default values doesn't match number of optional variables and inform the user
% (2) check we really need vargout here!?
% (3) default values doesn't have to be passed in a cell array they can be
% structure w/ same names 
% ------
% Code Info:
%   creation: 2014-11-07 by ShS -> shervin.safavi@gmail.com
%   modification: 
%       $ 201?

optionalVariableNames =  fieldnames(optionalVariables);

nInput = numel(optionalInputs);
nDefaultValues = numel(defaultValues);

for iInput = 1 : nInput
    if isempty(optionalInputs{iInput})
        optionalVariables.(optionalVariableNames{iInput}) = defaultValues{iInput};
    else
        optionalVariables.(optionalVariableNames{iInput}) = optionalInputs{iInput};
    end
end

for iVariable = nInput+1 : nDefaultValues
     optionalVariables.(optionalVariableNames{iVariable}) = defaultValues{iVariable};
end

varargout{1} = optionalVariables;




