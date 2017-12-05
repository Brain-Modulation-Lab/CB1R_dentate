function [center_overlap]=find_cencenoverlapped_combs(COPtS,COStP,comb)

tmp=intersect(COPtS{combs{s}(1)},COStP{combs{s}(1)});
for i=2:length(comb,2)
tmp=intersect(tmp,intersect(COPtS{comb(2)},COStP{comb(2)}));

end

tmp(ismember(tmp,ind))=[];
