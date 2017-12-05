function [output]=parse_stack(stack,crit_prc_ovl,class_confidence,zrange,SF,expos,Ccount)

% define the bins in the z-plane
nbins=40;
bins=0:nbins;
% normalization for z-stack size discrepancy
nfactor=max(bins)/size(stack.LipoMask,3);


ch_label=stack.ch_label;
ch_comb={stack.Stat(:).Combination};
[combs]=allchannelcombs(length(ch_label)-1);
% data normalization 
for i=1:length(ch_label)
    x=stack.Data.(ch_label{i});
    stack.Data.(ch_label{i})=x*expos(i)*SF(i);
    
end


for ch=1:length(ch_label)
    pri_mask=ch_label{ch};
    sec_masks=ch_label(setdiff(1:length(ch_label),ch));
    pri_ch=stack.Data.(pri_mask);
    
    sec_chs=cellfun(@(x) stack.Data.(x),sec_masks,'Uni',0);
    
    nsc=length(sec_masks);
    for s=1:nsc
        pri_sec(s)=find(contains(ch_comb,[pri_mask '__' sec_masks{s}]));
        sec_pri(s)=find(contains(ch_comb,[sec_masks{s} '__' pri_mask]));
    end
    
    % find objects with centroids within the selected z-range
    [~,~,whichbin]=histcounts(bsxfun(@times,stack.Stat(find(contains(ch_comb,pri_mask),1,'first')).CenZ,nfactor),bins);
    ind=find(ismember(whichbin,setdiff(1:nbins,zrange)));
    
    %%
    % find objects with percentage overlap
    i=1;
    POPtS=cell(length(sec_masks),1); COPtS=POPtS; COStP=COPtS;
    for c=pri_sec
        
        tmp=stack.Stat(c).PrcOverlap>crit_prc_ovl;
        [r,~]=ind2sub(size(tmp),find(tmp));
        POPtS{i,1}=unique(r);
        i=i+1;
    end
    %find objects with center-center overlap
    i=1;
    % objects overlapping the primary channel centroid
    for c=pri_sec
        tmp=stack.Stat(c).CenOverlap;
        [r,~]=ind2sub(size(tmp),find(tmp));
        COPtS{i,1}=r;
        i=i+1;
    end
    
    % objects overlapping the secondary channel centroid
    i=1;
    for c=sec_pri
        tmp=stack.Stat(c).CenOverlap;
        COStP{i,1}=vertcat(stack.Stat(c).OverlapObjID(logical(tmp)));
        i=i+1;
    end
    
    
    % find center center overlapped objects based on combinations
    CCOV=cell(length(combs),1);
    for s=1:length(combs)
        tmp=intersect(COPtS{combs{s}(1)},COStP{combs{s}(1)});
        for i=2:size(combs{s},2)
            tmp=intersect(tmp,intersect(COPtS{combs{s}(2)},COStP{combs{s}(2)}));
            
        end   
        % remove objects outside of the z-range
        tmp(ismember(tmp,ind))=[];        
        CCOV{s}=tmp;
    end
    % remove objects overlapping with objects from other channels
    for s=1:length(combs)
        sec_interct=setdiff(1:length(ch_label)-1,combs{s});
        for i=1:length(sec_interct)
            CCOV{s}=CCOV{s}(~ismember(CCOV{s},POPtS{sec_interct(i)}));
            
        end   
    end
    
    CCOVINT=cell(length(combs),1);
    for s=1:length(CCOV)
        if ~isempty(CCOV)
            CCOVINT{s}(:,1)= arrayfun(@(x) mean(pri_ch(stack.CC{ch}.PixelIdxList{x})),CCOV{s});
            for s2=1:length(sec_masks)
                CCOVINT{s}(:,s2+1)= arrayfun(@(x) mean(sec_chs{s2}(stack.CC{ch}.PixelIdxList{x})),CCOV{s});
                
            end
        end
        
    end
    [CCOV,CCOVINT]=classify_puncta(CCOV,CCOVINT,stack,POPtS,pri_ch,sec_chs,ch,class_confidence);
    
end

%%
% output(1,:)=[{'comb'}  {'density_per_volume'} {'density_per_cell'} {'count'} {'size'} {pri_mask} sec_masks(:)' {'puncta'}];
% for s=1:length(CCOV)
%     tmp=[CCOV{s} bsxfun(@times,CCOVINTS{s},expos(1:4).*SF(1:4))];
%     x=[];
%     for y=1:length(combs{s})
%         x=cat(1,x,sec_masks{combs{s}(y)});
%     end
%     output{s+1,1}=x;
%     if ~isempty(CCOV{s})
%
%         output(s+1,6:9)=num2cell(mean(tmp(:,2:end)));
%         output(s+1,2)=num2cell(length(CCOV{s})/vol);
%         output(s+1,3)=num2cell(length(CCOV{s})/Ccount);
%         output(s+1,4)=num2cell(length(CCOV{s}));
%
%         switch s
%             case {1 5 7}
%                 id=CCOV{s};
%                 tmp2=zeros(length(id),1);
%                 for i=1:length(id)
%                     tmp2(i)= max(Stat(4).vol(nonzeros(Stat(1).OverlapObjID(id(i),:))));
%                 end
%                 output(s+1,5)=num2cell(median(tmp2)*(.267*0.267*0.25));
%
%             otherwise
%
%                 output(s+1,5)=num2cell(median(Stat(1).vol(CCOV{s}))*(.267*0.267*0.25));
%         end
%         output{s+1,10}=tmp;
%     else
%         output(s+1,6:9)=num2cell(NaN(1,4));
%         output(s+1,2)=num2cell(0);
%         output(s+1,3)=num2cell(0);
%         output(s+1,4)=num2cell(0);
%         output(s+1,5)=num2cell(NaN);
%         output{s+1,10}={0};
%
%     end
% end


end

