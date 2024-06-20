% Open file selection dialog and allow selection of multiple data files
cd('data')
[fileNames, filePath] = uigetfile({'*.csv'; '*.xlsx'}, 'MultiSelect','on');

% Check if any files were selected
if iscell(fileNames)
    numFiles = length(fileNames);
else
    numFiles = 1;
    fileNames = {fileNames};
end

% initialize a structure to store the results
results = struct();

% Process each selected file
for i = 1:numFiles
    fullFilePath = fullfile(filePath, fileNames{i});
    data = readmatrix(fullFilePath);
    
    numSubjects = size(data,1); % sample size
    numCh = size(data,2);   % number of channels
    
    meffMean = zeros(numSubjects, 1);
    meffStd = zeros(numSubjects, 1);
    meffResample = zeros(numSubjects, 1000);
    
    % Calculate Meff 1000 times each N
    for subjectsIdx = 3:numSubjects
        for tmp = 1:1000
	        sampledData = randperm(numSubjects,subjectsIdx);
            Data = data(sampledData,:);
            realEig=real(eig(corrcoef(Data)));
            dEig=sort(realEig,'descend');
            meffResample(subjectsIdx, tmp)=real(((sum(sqrt(dEig)))^2)/sum(dEig));
        end
        meffMean(subjectsIdx)=mean(meffResample(subjectsIdx,:));
        meffStd(subjectsIdx)=std(meffResample(subjectsIdx,:));
    end

    % Store the results in the structure
    results(i).fileName = fileNames{i};
    results(i).meffMean = meffMean;
    results(i).meffStd = meffStd;
end


