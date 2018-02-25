% find all files to read 
files=dir('*.mat');
files=vertcat({files.name});
% parse out subjects 
sloc=strfind(files,'_');
sloc=vertcat(sloc{:});
sloc=sloc(:,1)';
subjects=unique(arrayfun(@(x,y) x{:}(1:y-1),files,sloc,'Uni',0));
%% define setting for the analysis 
round=1;
% you need to load the corresponding imaging_parameters_round and
% cell_counts mat files before running the analysis
% z-planes to be used 
zrange=5:35;
% threshold for mask overlap percentage to be considered real 
crit_prc_ovl=0;
% threshold for classification confidence to be considered real 
class_confidence=0.20;

% stack volume in microns
        vol=0.25*(zrange(2)-zrange(1))*(512*0.265)^2;
for sn=1:length(subjects)
    % find files for each subjects 
    sub_files=files(contains(files,[subjects{sn} '_']));
    fn=[];
    tic
    for f=1:length(sub_files)
        
        fname=sub_files{f};
        %load in data 
        stack = load(fname);
        % import the scaling factor and exposure time for that stack to
        % normalize the intensities, make sure you load the appropriate
        % imaging_parameters matfile for that round
        uloc=strfind(fname,'_');
        idx=contains(cellstr(StackName),fname(1:uloc(2)-1));
        expos=ExposureTimesNorm(idx,:);
        
        SF=ScalingFactorsNorm(idx,:);
        
       idx=contains(cellstr(Site_name),fname(1:uloc(2)-1));
    if sum(idx)>0
        Ccount=CellCounts(idx);
    else 
    Ccount=NaN;    
        
    end 


        
        [output]=parse_stack(stack,crit_prc_ovl,class_confidence,zrange,SF,expos,Ccount,vol);
    
       lipoInt=bsxfun(@times,stack.LipoMask,stack.Data.Lipofuscin)*expos(5)*SF(5);
        
        SiteStat(f).site=fname;
        SiteStat(f).puncta=output;
        SiteStat(f).count=Ccount;
        SiteStat(f).vol=vol;
        SiteStat(f).Lipo=[sum(lipoInt(:)) sum(stack.LipoMask(:))];
        
        
        
        
    end
    toc
    FinalStat(sn).Subject=subjects{sn};
    FinalStat(sn).Results=SiteStat;
    clearvars SiteStat
    
end
