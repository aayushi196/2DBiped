function [rmse,actY] = calculateRMSE(fis,x,y)

% Specify options for FIS evaluation
persistent evalOptions
if isempty(evalOptions)
    evalOptions = evalfisOptions("EmptyOutputFuzzySetMessage","none", ...
        "NoRuleFiredMessage","none","OutOfRangeInputValueMessage","none");
end

% Evaluate FIS
actY = evalfis(fis,x,evalOptions);

% Calculate RMSE 
del = actY - y;
rmse = sqrt(mean(del.^2));

end