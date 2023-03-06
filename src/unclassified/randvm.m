function R = randvm(kappa, varargin) 
% R = randvm(kappa, dim, mu) 
% This function returns pseudorandom values drawn from  a von Mises
% distribution. User has to specify the concentration parameter "kappa" (a
% a reciprocal measure of dispersion). It's important to note that 1/kappa
% is analogous to sigma^2 for a normal distribution and for kappa = 0, 
% distribution is uniform, for large kappa all the samples will be
% concentrate around the mean (mu: default is zero). 
% 
% For more details see: 
% [1] Fisher NI, Statistical analysis of circular data (2000); sec. 3.3.6
% [2] https://en.wikipedia.org/wiki/Von_Mises_distribution
% [3] http://mathworld.wolfram.com/vonMisesDistribution.html
%
% EXAMPLE:
% an example syntax (if it's needed)
%
% ------
% Input:
% 1     kappa: scalar
%       expnation
% (2)   dim: scalar or a vector
%       expnation
% (3)   mu: scalar 
%       (default value = 0)
%       an angle between 0 and 2*pi
% 
% Output:
% 1     R: 
%
% ------
% potential improvments:
% (1) if the required dimention is more than 2 the code did not incorporate
% it
% ------
% Code Info:
%   creation: 2018-04-08 by SS (shervin.safavi@tuebingen.mpg.de)
%   modification:
%       $ 2018-??-?? ???
% ------
% see also randn

%% Handle optional inputs (varargin):
optionalVariables.dim   = [];     defaultValues{1} = 1;
optionalVariables.mu	= [];     defaultValues{2} = 0;

optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

% decide about the dimntion of the output R
if numel(optionalVariables.dim) == 1    
    R_nonWrapped = nan(optionalVariables.dim, optionalVariables.dim);
    % (trying to be consistent with other MATLAB pseudorandom generator routins)
else
    R_nonWrapped = nan(optionalVariables.dim(1), optionalVariables.dim(2));
end
    
%% pseudorandom angle generation
% variable name are chosen according to section 3.3.6 (page 49) of Fisher
% 2000 [1]

a = 1 + (1 + 4 * kappa ^ 2) ^ .5;
b = (a - (2 * a) ^ .5) / (2 * kappa);
r = (1 + b ^ 2) / (2 * b);

% following the iterative algorithm explained in section 3.3.6 (page 49)
% Fisher 2000 [1]:

nR = numel(R_nonWrapped);       % number random angles
for iR = 1 : nR
  while true
      u = rand(1,3); 
      % U(1), U(2) and U(3) correspond to U_1, U_2 and U_3 in [1]
      % U_1, U_2 and U_3 are random numbers drawn from a uniform
      % distribution between 0 and 1

      z = cos(pi * u(1));                   % related to step 1
      f = (1 + r * z) / (r + z);            % related to step 1
      c = kappa * (r - f);                  % related to step 1

      if c * (2 - c) - u(2) > 0             % related to step 2                      
        break
      elseif log(c / u(2)) + 1 - c >= 0     % related to step 3
        break
      end      
  end
  R_nonWrapped(iR) = sign(u(3) - 0.5) * acos(f) + optionalVariables.mu;  
                                            % related to step 4
end

R = angle(exp(1i * R_nonWrapped));          % related to step 4


