% fis = readfis('E:\Work\projects\biped\fuzzy_logic_library\Navigation-With-Fuzzy-Logic-master\Fuzzy\Triangular.fis');
% fis = readfis('Gaussian.fis');
% fis
% fis.Inputs(1)
% fis.Inputs(2)
% fis.Inputs(1).MembershipFunctions
% fis.Inputs(2).MembershipFunctions
% fis.Rules
% plotfis(fis)
% plotmf(fis,'input',1)
% plotmf(fis,'input',2)
% plotmf(fis,'output',1)
% gensurf(fis)
% evalfis(fis,[4 3])

X(:,1) = err_1;
X(:,2) = err_dot_1;
Y = corr_1;
trnX = X(1:2:end,:); % Training input data set
trnY = Y(1:2:end,:); % Training output data set
vldX = X(2:2:end,:); % Validation input data set
vldY = Y(2:2:end,:); % Validation output data set
dataRange_input = [min(X)' max(X)'];
dataRange_output = [min(Y)' max(Y)'];
fisin = mamfis;
name = ['error', 'derror','correction'];
for i = 1:2
    fisin = addInput(fisin,dataRange_input(i,:),'Name',name(i),'NumMFs',7);
end
fisin = addOutput(fisin,dataRange_output(1,:),'Name',name(3),'NumMFs',7);
figure
plotfis(fisin)

options = tunefisOptions('Method','particleswarm',...
    'OptimizationType','learning', ...
    'NumMaxRules',50);
options.MethodOptions.MaxIterations = 20;
rng('default')
runtunefis = true;
if runtunefis
    [fisout1 optimout1]= tunefis(fisin,[],trnX,trnY,options); 
else
    tunedfis = load('tunedfismpgprediction.mat'); %#ok<UNRCH>
    [fisout1 optimout1] = tunedfis.fisout1;
    fprintf('Training RMSE = %.3f MPG\n',calculateRMSE(fisout1,trnX,trnY));
end
plotfis(fisout1)
[fisout1.Rules(1:3).Description]' 
plotActualAndExpectedResultsWithRMSE(fisout1,vldX,vldY)

[in,out,rule] = getTunableSettings(fisout1);
options.OptimizationType = 'tuning';
options.Method = 'patternsearch';
options.MethodOptions.MaxIterations = 60;
if runtunefis
    rng('default')
    fisout = tunefis(fisout1,[in;out;rule],trnX,trnY,options);
else
    fisout = tunedfis.fisout; %#ok<UNRCH>
    fprintf('Training RMSE = %.3f MPG\n',calculateRMSE(fisout,trnX,trnY));
end
figure
plotfis(fisout)
plotActualAndExpectedResultsWithRMSE(fisout,vldX,vldY);