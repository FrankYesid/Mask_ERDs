load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
set(0,'DefaultFigureWindowStyle','docked')
for fold = 1
    for cl = 1:2
        for ch = [8,10,12,20]
            dat = erd{fold}{cl}.ERDS{ch}.erds;
            figure
            imagesc((1:1751)./fs,1:17,dat',[-1 1.5])
            yticks(1:size(filter_bank,1));
            yticklabels({num2str(filter_bank)})
            colormap(erdcolormap)
            colorbar
            axis xy;
            title(['Sub: ' num2str(8) ' Cl: ' num2str(cl) ' Ch: ' num2str(ch)])
        end
    end
end