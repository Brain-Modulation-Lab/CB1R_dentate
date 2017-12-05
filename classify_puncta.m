function [CCOV,CCOVINT]=classify_puncta(CCOV,CCOVINT,stack,POPtS,pri_ch,sec_chs,ch,class_confidence)
        class=[];
        for i=1:length(CCOV)
            class=cat(1,class,i*ones(length(CCOV{i}),1));
        end
        training_data=vertcat(CCOVINT{:});
        if ~isempty(class)
            % build the naive bayesian classifier
            mdl=  fitcnb(training_data, class,'Prior','Uniform','Distribution','kernel');
        end
 
        if ~isempty(mdl)
            % get unclassified puncta
            test=unique(vertcat(POPtS{:}));
            test(ismember(test,vertcat(CCOV{:})))=[];
            test_int=zeros(length(test),length(ch_label));
            test_int(:,1)= arrayfun(@(x) mean(pri_ch(stack.CC{ch}.PixelIdxList{x})),test);
            for s2=1:length(sec_masks)
                test_int(:,s2+1)= arrayfun(@(x) mean(sec_chs{s2}(stack.CC{ch}.PixelIdxList{x})),test);
                
            end
            
            [cpre,~,Cost] = predict(mdl,test_int);
            for co=1:size(Cost,2)
                test_class=cpre.*bsxfun(@le,Cost(:,co),class_confidence);
                if ~isempty(nonzeros(test_class))
                    CCOV{unique(nonzeros(test_class))}=cat(1,CCOV{unique(nonzeros(test_class))},test(bsxfun(@gt,test_class,0)));
                    CCOVINT{unique(nonzeros(test_class))}=cat(1,CCOVINT{unique(nonzeros(test_class))},test_int(bsxfun(@gt,test_class,0),:));
                    
                end
            end
        
        
    end