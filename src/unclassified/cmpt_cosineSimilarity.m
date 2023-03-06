function cosSim = cmpt_cosineSimilarity(v1,v2)
% cosSim = cmpt_cosineSimilarity.m(v1,v2)

% check the input
if isvector(v1) == 0 || isvector(v2) == 0
    error('Input to the function should be vectors')
end
if length(v1)~=length(v2)
    error('v1 and v2 need to have identical length!')
end

% compute the cosine similarity
dotProduct   = dot(v1, v2);
v1_norm   = norm(v1);
v2_norm   = norm(v2);
% normProduct = nx*ny;
if dotProduct == 0
    cosSim = 0;
else
    cosSim   = dotProduct / ((v1_norm) * (v2_norm));
end