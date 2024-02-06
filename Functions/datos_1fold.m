clc
for s = [1:28,30:33,35:52]
    fprintf(['Subject: ',num2str(s),'\n'])
    load(['I:\ERD_giga\MI\30_folds_training\ERD_folds_s',num2str(s),'.mat'])
    erd_o = cell(1,30);
    erd_ = cell(1,30);
    for fold = 1:30
        tem = ERDsfilt_o{fold};
        ma_ = max(cell2mat(cellfun(@(y) max(cell2mat(cellfun(@(x) max(x(:)),y,'UniformOutput',false))),tem,'UniformOutput',false)));
        mi_ = min(cell2mat(cellfun(@(y) min(cell2mat(cellfun(@(x) min(x(:)),y,'UniformOutput',false))),tem,'UniformOutput',false)));
        erd_o{fold} = cellfun(@(y) cellfun(@(x) 1.5+(((x-mi_).*(-1-(1.5)))./(ma_-mi_)),y,'UniformOutput',false),tem,'UniformOutput',false);
        clear tem
        tem = ERDsfilt_{fold};
        ma_ = max(cell2mat(cellfun(@(y) max(cell2mat(cellfun(@(x) max(x(:)),y,'UniformOutput',false))),tem,'UniformOutput',false)));
        mi_ = min(cell2mat(cellfun(@(y) min(cell2mat(cellfun(@(x) min(x(:)),y,'UniformOutput',false))),tem,'UniformOutput',false)));
        for cl = 1:2
            for ch = 1:64
                te = tem{cl}{ch};
                te = 1.5+(((te-mi_).*(-1-(1.5)))./(ma_-mi_));
                ta = tem{cl}{ch};
                te(ta==0)=0;
                erd_{fold}{cl}{ch} = te;
            end
        end
    end
    clear ERDsfilt_ ERDsfilt_o
    ERDsfilt_ = erd_;
    ERDsfilt_o = erd_o;
    clear erd_ erd_o
    save(['D:\ERD_folds_s',num2str(s)],'ERDsfilt_','ERDsfilt_o')
    clear ERDsfilt_ ERDsfilt_o
end
clc