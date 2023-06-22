function [Pi_f,mask_i] = fncmask(rho,threshold,filter_bank,s,X,graf,sub)
col = [215,215,215]; col2 = [014,041,075]; %oscuro
c = [linspace(col2(1),col(1),64)',linspace(col2(2),col(2),64)',...
    linspace(col2(3),col(3),64)']/255;

if sum(isnan(rho))>=1
    rho(isnan(rho)==1) = threshold(end);
end %para los nan detectados
mask_i = reshape(rho,[size(filter_bank,1) size(X{1},2)]);

% promedio correlaciones by channel
Pi_c = zeros(1,size(X{1},2));
for i=1:22; Pi_c(i) = mean(mask_i(:,i));end

% promedio correlaciones by frequency
Pi_f = zeros(1,size(filter_bank,1));
for i=1:17; Pi_f(i) = mean(mask_i(i,:)); end

%% grafica
if graf == 1
    subplot(2,5,sub)
%     figure; 
    x = 0.25;     y1 = 0.15;     w = 0.55;      h = 0.8;
    %     subplot('position',[x,y1,w,h]);
    imagesc(mask_i,[0,1]); %min(rho),max(rho)
%     title(num2str(s),'FontSize',15,'Interpreter','latex')
    title(['Fold:',num2str(sub)],'FontSize',15,'Interpreter','latex')
    %         colormap(c)
    axis xy;  xticks([]);
    xcol = x+w+0.003+0.01+0.06+0.01;
%     colorbar('Position',[xcol y1 0.01 h]);
    
    if sub == 1 || sub == 6
    yticks(1:size(filter_bank,1));
    yticklabels({num2str(filter_bank)})
    ylabel('Frequency','FontSize',36,'Interpreter',  'latex')
    else
         yticklabels({})
    end
    
%     x = 0.25; y2 = 0.12;  w = 0.55; h1 = 0.10;
%     subplot('position',[x,y2,w,h1])
%     bar(1:22,Pi_c,'FaceColor',c(end,:),...
%         'EdgeColor',c(end,:)); xlim([0.6,22.4]);ylim([0,1])
    if sub == 6 || sub == 7 || sub == 8 || sub == 9 || sub == 10
    xticks(1:size(X{1},2));
    xticklabels({1:7,'C3',9,'Cz',11,'C4',13:22})
    xlabel('Channels','FontSize',36,'Interpreter','latex')
    end
    % marginal en el eje y
%     x_m2 = x+w+0.003+0.01; y_3 = y1;%%%%
%     w_m2 = 0.06; h_m2 = h;
%     subplot('position',[x_m2,y_3,w_m2,h_m2])
%     barh(1:17,Pi_f,'FaceColor',c(end,:),'EdgeColor',c(end,:));
%     ylim([0.6,17.4]);xlim([0,1])
%     yticks([])
        saveas(gcf,[SUBJECTS_DIR2 filesep 'masks' filesep SUBJECTS{s}],'epsc')
            saveas(gcf,[SUBJECTS_DIR2 filesep 'masks' filesep SUBJECTS{s}],'fig')
end
end