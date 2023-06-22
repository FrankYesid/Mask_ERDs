dir = 'F:\BCI\BCICIV_2a_08\ERD_';
erd_ = zeros(22,1751,17);
 
% set(0,'DefaultFigureWindowStyle','docked')
for s = 8
    load(['F:\New_ERD_abri_2019\ERD_folds_s' num2str(s) '.mat'])
    for clas = 1
        for fold = 1
            load([dir num2str(s) '_' num2str(clas) '_' num2str(fold) '_.mat'])
            r = ERD;
            for chn = 1:length(r.ERDS)
                sig_matrix = (r.ERDS{chn}.cl > 0 & r.ERDS{chn}.cu > 0) | ...
                    (r.ERDS{chn}.cl < 0 & r.ERDS{chn}.cu < 0);
                r.ERDS{chn}.erds = sig_matrix .* r.ERDS{chn}.erds;
                ERD_p(chn,:,:) = r.ERDS{chn}.erds;
                ERD_m(chn,:,:) = erd{fold}{clas}.ERDS{chn}.erds;
            end;            
        end
    end
    
end

load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
figure;
imagesc((1:1751)/250,1:17,squeeze(ERD_m(8,:,:))',[-1 1.5])
axis xy; colormap(erdcolormap); colorbar
figure
imagesc((1:1751)/250,1:17,squeeze(ERD_p(8,:,:))',[-1 1.5])
axis xy; colormap(erdcolormap); colorbar

