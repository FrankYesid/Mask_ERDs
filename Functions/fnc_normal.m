clear; close all; clc
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

x = [-3:.1:3];
norm = normpdf(x,0,1);

plot(x,norm)

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'XTickLabel','','YTickLabel','');
% Set the remaining axes properties
set(axes1,'FontSize',5,'XColor',[0 0 0],'XTickLabel','','YColor',[0 0 0],...
    'YTickLabel','','ZColor',[0 0 0]);
% Create line
annotation(figure1,'line',[0.51890756302521 0.51890756302521],...
    [0.919275123558484 0.113673805601318]);

% Create line
annotation(figure1,'line',[0.882352941176441 0.882352941176441],...
    [0.126853377265239 0.112026359143328]);

% Create line
annotation(figure1,'line',[0.152310924369748 0.152310924369748],...
    [0.125205930807249 0.110378912685338]);

% Create textbox
annotation(figure1,'textbox',...
    [0.871798319327729 0.070840197693575 0.0284117647058845 0.0428336079077429],...
    'String','\sigma_{l}',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.137554621848738 0.0724876441515651 0.0231596638655475 0.042833607907743],...
    'String','\sigma_{u}',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.507302521008403 0.0675453047775948 0.0368151260504203 0.042833607907743],...
    'String','\mu',...
    'FitBoxToText','off',...
    'EdgeColor','none');

%%
saveas(gcf,'Figure','epsc')