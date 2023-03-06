function kappa = MLest_kappa_forVonMisesBasedPLV1(plv, n)
% similar/exact to "estimate_kappa_forVonMisesBasedPLV1" in cribay
% kappa = estimate_kappa_forVonMisesBasedPLV1(plv)

if n < 15
    kappaEst_ML = MLestimate_kappa(plv);
    if kappaEst_ML < 2
        kappa = max([kappaEst_ML 0]);
    else
        kappa = kappaEst_ML * (n-1) / (n^3 + n);
    end 
        
else % if n >= 15
    kappa = MLestimate_kappa(plv);
    
end

end

function kappaEst_ML = MLestimate_kappa(plv)
if plv < 0.53
    kappaEst_ML = 2*plv + plv^3 + 5*plv^5/6;
    
elseif plv >= 0.85
    kappaEst_ML = 1 / (plv^3 - 4*plv^2 + 3*plv);

else
    kappaEst_ML = -.4 + 1.39*plv + .43/(1 - plv);

end
end