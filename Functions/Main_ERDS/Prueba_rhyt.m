clear ; clc
band  ='beta';
for s = 8
    for len = 2%[2,1.5,1,0.5,0.2]
        for Distan = 4%[1,3,4]
            for fold = 1%:Nfolds
                Dis1 = load(['C:\Users\Lufe\Desktop\Prueba_dif',filesep,'2Sub_',num2str(s),'_w',num2str(len*1000),'_Met',num2str(Distan),'_fold',num2str(fold),band]);
                Dis2 = load(['C:\Users\Lufe\Desktop\Prueba_dif',filesep,'2Sub_',num2str(s),'_w',num2str(len*1000),'_Met',num2str(Distan),'_fold',num2str(fold),band,'2']);
                nchan = size(Dis1.Disaw_{1},1);
                nw = size(Dis1.Disaw_{1},2);
                for cl = 1:2
                    switch Distan
                        case 1
                            %% Method 1: inner product
                            %         disp('Method 1: inner product')
                            %         Dis_ = zeros(nchan,nfreq,nw);
                            for ch = 1:nchan
                                for tao = 1:nw % abs(Disa_{cl}(:,tao))/max(abs(Disa_{cl}(:,tao)))
                                    Dis_(ch,tao) = abs(Dis1.Disaw_{cl}{ch,tao}*Dis2.Disaw_{cl}{ch,tao}'); %norm(ERD_{folds}{1}{ch}{fre}{tao}-ERD_{folds}{2}{ch}{fre}{tao});
                                end
                            end
                            % Normalization
                            % Dis_ = Dis_./max(max(max(Dis_)));
                            
                        case 3
                            %% Method 3: Euclidean
                            %         disp('Method 3: Euclidean')
                            %         Dis_ = zeros(nchan,nfreq,nw);
                            for ch = 1:nchan
                                for tao = 1:nw
                                    Dis_(ch,tao) = norm(Dis1.Disaw_{cl}{ch,tao}-Dis2.Disaw_{cl}{ch,tao},2);
                                end
                            end
                            
                            
                        case 4
                            %% Method 4: Correlation
                            %         disp('Method 4: Correlation')
                            for ch = 1:nchan
                                for tao = 1:nw
                                    Dis_(ch,tao) = abs(corr(Dis1.Disaw_{cl}{ch,tao}',Dis2.Disaw_{cl}{ch,tao}'));
                                end
                            end
                            
                        otherwise
                            disp('Distance 1) inner product, 2) 1-inner product, 3) Euclidean, 4) Correlation, 5) 1-Correlation')
                    end
                    Disa_{cl} = Dis_./max(max(max(Dis_)));
                end
                %                 save(['C:\Users\Lufe\Desktop\Prueba_dif',filesep,'Dis_2Sub_',num2str(s),'_w',num2str(len*1000),'_Met',num2str(Distan),'_fold',num2str(fold),band],'Disa_')
            end
        end
    end
end


%  Plot de las figuras
% clear;clc;
band  ='beta';
load('G:\Dropbox\ERD\Codes\Topoplots\BCICIV_2a\electrodesBCICIV2a.mat')
%load('BCICIV_2a\labels.mat')
%load('BCICIV_2a\layout.mat')
load('G:\Dropbox\ERD\Codes\Topoplots\HeadModel.mat')   % model of the head.
Fil = 3; Col = 3;              % topoplot en subplots para cada uno de los sujetos.
t_e = 70;                        % tamaño visual de los electrodos seleccionados.
M1.xy = elec_pos(:,1:2);% posicion de los canales.
M1.lab = Channels;        % nombre de los canales.
% mask = reshape(rho,[17 22]) <= thresh; %
sel = 1:22;                     %find(mask(band,:)==1); % Orden de los canales segun una relevancia.
% relevancia para los canales la variable es w.
% tmp = W{fold,band}(caract,:); % rel = linspace(1,10,numel(M1.lab));     % para colocar mas importancia al electrodo.
% rel = ones(1,22);

% for sss = 8
%     for cl = 1%:2
%         for len = 2 %[2,1.5,1,0.5,0.2]
%             for Distan = [4]
%                 for fold = 1
%                     load(['C:\Users\Lufe\Desktop\Prueba_dif',filesep,'Dis_2Sub_',num2str(sss),'_w',num2str(len*1000),'_Met',num2str(Distan),'_fold',num2str(fold),band])
%                     rel = Disa_{cl};
%                     fs = 250;
%                     nsamples = 1751;
%                     over = 0.1;
%                     t = 1:round(len*fs*over):(nsamples-(len*fs));
%                     imagesc(t./250,1:22,rel,[0 1])
%                     axis xy
%                     xlim([0 7])
%                 end
%             end
%         end
%     end
% end

for sss = 8
    for cl = 1%:2
        for len = 2 %[2,1.5,1,0.5,0.2]
            for Distan = [4]
                for fold = 1
                    %                     load(['C:\Users\Lufe\Desktop\Prueba_dif',filesep,'Dis_2Sub_',num2str(sss),'_w',num2str(len*1000),'_Met',num2str(Distan),'_fold',num2str(fold),band])
                    for tao = [1,5,8,10,13,14,18,22,25]
                        rel = abs(Disa_{cl}(:,tao));%/max(abs(Disa_{cl}(:,tao)));
                        %                         rel = zeros(1,22);
                        figure;
                        fig = gca;
                        %                 fig.Parent.PaperPosition = [2.91 4.15 2.69 2.69];
                        %                         fig.Parent.OuterPosition = [1998 541 406 468];
                        fig.Parent.OuterPosition = [2311 678 243 315];
                        %                         rel = squeeze(er{cl}(:,t_(tim),fre));
                        pos = M1.xy;
                        label = M1.lab;
                        % tam - tamaño de la cabeza y posición de los electrodos.
                        tam = 1; tam2 = 0.96; tam3 = 0.98;
                        for i=1:2
                            pos(:,i) = tam2.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
                        end
                        xc = HeadModel(1,:);
                        yc = HeadModel(2,:);
                        %                         load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
                        hold on
                        % Topoplot
                        x = pos(:,1);
                        y = pos(:,2);
                        tmp = [x,y,x*0]*rotz(2);
                        x = tmp(:,1);
                        y = tmp(:,2);
                        pos(:,1) = x;
                        pos(:,2) = y;
                        GS = 500;
                        xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
                        yi = linspace(min(y)-0.05, max(y)+0.05, GS);       % y-axis for interpolation (row vector)
                        [Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
                        %% Creating data mask
                        [TH R] = cart2pol(Xi,Yi);
                        Zi(R>0.5) = NaN;
                        deltax = xi(2)-xi(1); % length of grid entry
                        deltay = yi(2)-yi(1); % length of grid entry
                        max(rel)
                        h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
                        %                         shading interp % interpola colores.
                        %                 scatter(tam3.*x(sel),(tam3+0.04).*y(sel),t_e,'b','filled')
                        %                 text(pos(sel,1)*tam3-0.02,pos(sel,2)*tam3+0.02,label(sel),'Interpreter',...
                        %                     'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        %                     'FontSize',10);
                        plot(tam*xc,tam*yc+0.001,'k','LineWidth',3)
                        %                         colormap(erdcolormap)
                                                caxis([0 1])
                        axis off
                        hold off
                        axis square
                        axis image
                        fs = 250;
                        nsamples = 1751;
                        over = 0.1;
                        t = 1:round(len*fs*over):(nsamples-(len*fs));
                        title([num2str(t(tao)/250) '-' num2str((t(tao)/250+len))])
                        %                 fig = gca;
                        %                 saveas(fig,['D:\Sub_beta_up',num2str(sss),'_c_',num2str(cl),'_t_',num2str(tim),'_f_',num2str(fre),'__'],'png')
                        %                         print('-depsc2', ['I:\Sub_',num2str(sss),'_c_',num2str(cl),'_t_',num2str(tim),'_f_',num2str(fre),'__.eps'])
                        %                 view([-90 90])
                        %                         close
                    end
                end
            end
        end
    end
end

