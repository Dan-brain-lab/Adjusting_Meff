clear;clc
cd('data');

% Define target values for different tasks, fill in with actual values
% targetValues.Stroop = ;
% targetValues.WTT = ;
% targetValues.GoNogo = ;
% targetValues.Placebo = ;

targetValues.IQ_low_deoxy_ch = 24.70136085;


% File selection interface
[fileName, filePath] = uigetfile({'*.csv'; '*.xlsx'});

% Exit if no file is selected 
if isequal(fileName, 0)
    disp('User pressed cancel');
    return;
end

fullFilePath = fullfile(filePath, fileName);
data = readmatrix(fullFilePath);

% Retrieve target value for the task
[~, task, ~] = fileparts(fileName);

if isfield(targetValues, task)
    targetValue = targetValues.(task);
else
    warning('Actual value for %s is not defined. Using default value.', task);
    targetValue = 0; 
end

[numSubjects, numCh] = size(data);
sampleRange = 10:numCh;
iterations = 1000;

% Define the fit type and option
ft = fittype('a*exp(b*x) + c', 'independent', 'x');
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.MaxFunEvals = 100000;
opts.StartPoint = [-20 0 20];

results = struct(); 

% Loop over each N to perform resampling and fitting
for N = sampleRange
    x = 3:N; %Define the range of data points used for fitting
    meffResample = zeros(length(x), iterations);
    
    params = zeros(iterations,3);
    predictedMeffs = zeros(iterations,1);

    for iter = 1:iterations
        isValid = false;

        % Resample the data and calculate Meff for each N
        while ~isValid
            for i = 3:N
                sampleData = randperm(numSubjects,i);
                realEig = real(eig(corrcoef(data(sampleData,:))));
                dEig = sort(realEig, 'descend');
                meffResample(i-2,iter) = real(((sum(sqrt(dEig)))^2)/sum(dEig));
            end

            if any(isnan(meffResample(:, iter)))
                continue;
            else
                isValid = true;
            end
        end

        % Fit the model to the resampled Meff
        [fitResult, gof] = fit(x', meffResample(:,iter), ft, opts);
        predictedMeff = fitResult.a * exp(fitResult.b * (numCh + 1)) + fitResult.c;
        predictedMeffs(iter,1) = predictedMeff;
        params(iter, :) = [fitResult.a, fitResult.b, fitResult.c];
        
    end

    % Store the results for each N
    results(N).subjects = N;
    results(N).predictedMeffsMean = mean(predictedMeffs);
    results(N).predictedMeffsStd = std(predictedMeffs);
    results(N).meanDiff =  mean(abs(predictedMeffs - targetValue));
    results(N).stdDiff = std(abs(predictedMeffs - targetValue));
end

% Change directory to create the results file
cd('../results');
fields = fieldnames(results);

filename = ['results_' task '.csv']; 
fid = fopen(filename, 'w');
fprintf(fid, 'N, Average of predicted Meff, SD of predicted Meff, Average of difference, SD of differece\n');
for i = sampleRange
    fprintf(fid, '%d,', i);
    fprintf(fid, '%f,', results(i).predictedMeffsMean);
    fprintf(fid, '%f,', results(i).predictedMeffsStd);
    fprintf(fid, '%f,', results(i).meanDiff);
    fprintf(fid, '%f\n', results(i).stdDiff);
end
fclose(fid);
