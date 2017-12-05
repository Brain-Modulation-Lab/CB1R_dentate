%% Import data from spreadsheet
files=dir('*.xlsx');

%% Import the data
StackName=[];  ScalingFactors=[];  ExposureTimes=[];
for i=1:length(files)
[~, ~, data] = xlsread(files(i).name,'Sheet1');

data(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),data)) = {''};

idx = cellfun(@ischar, data);
data(idx) = cellfun(@(x) string(x), data(idx), 'UniformOutput', false);
StackName=cat(1,StackName,data(1:5:end,1));  
tmp=cellstr(data(:,2));
idx1=strfind(tmp,': ');
idx2=strfind(tmp,']');
SF=cell2mat(cellfun(@(x,y,z) str2double(z(x+1:y-1)),idx1,idx2,tmp,'Uni',0));
ScalingFactors=cat(1,ScalingFactors,reshape(SF,5,[])'); 
ExposureTimes=cat(1,ExposureTimes,reshape(cell2mat(data(:,3)),5,[])');

end
%% Clear temporary variables
ScalingFactorsNorm=ScalingFactors./min(ExposureTimes);
ExposureTimesNorm=ExposureTimes./min(ExposureTimes);
save('imaging_parameters_round1.mat','StackName','ScalingFactors','ExposureTimes','ScalingFactorsNorm','ExposureTimesNorm')