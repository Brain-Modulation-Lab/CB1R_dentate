function [combs]=allchannelcombs(ch)
% create a cell array of all possible channel combination 
combs=[];
tmp=perms(1:ch);
for i=1:ch
combs=cat(1,combs,num2cell(unique(sort(tmp(:,1:i),2,'ascend'),'rows'),2));  
 
end

end