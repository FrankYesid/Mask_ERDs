clear all; clc % limpiar datos.

%% load data
load('BCICIV_2a.mat');

% ERD con Wavelet band
t = 0:1/250:7-(1/250);
t1 = 0;
t2 = 2;
lv = 5;
subjeto = 9;
name = 'adwd3';
Y_suj = cell(9,1);

for su = subjeto%1:9
    Sujeto = X{su};
    etiqueta = y{su};
    
    % escojer los train de las etiquetas 1 y 2
    Data = Sujeto(ismember(etiqueta,[2]));
   
    [WaveRec  WaveCoef]= MultiWPdec(Data,'sym15',lv);
    [BestRec,BestT] = f_BestTree(WaveRec,WaveCoef);
    Y_suj{su}(:,:,:,:) = BestT;  % (scale,channel,time,trials)
    %     end
end

%% tiempo de referencia
temp1 = abs(t - t1);
min1 = min(temp1);
temp2 = abs(t - t2);
min2 = min(temp2);
ul = find(temp1 == min1);
up = find(temp2 == min2);

%
fprintf(['Et{Pnc(t,f):te[Ta,Tb]}\n'])
for su = subjeto
    fprintf(['sujeto: ' num2str(su) ' de ' '9' '\n'])
    % Et{Pnc(t,f):te[Ta,Tb]}
    r_nc{su,1} = squeeze(mean(Y_suj{su}(:,:,ul:up,:),3));
end

%
fprintf(['En{r_nc(f)}\n'])
for y = subjeto
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % En{r_nc(f)}
    r_c{y,1} = squeeze(mean(r_nc{y,1},3));
end
%
fprintf(['En{Pnc(t,f)}\n'])
for y = subjeto
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % En{Pnc(t,f)}
    m_c{y,1} = squeeze(mean(Y_suj{y},4));
end

fprintf(['ERD = m_c/r_c\n'])
for y = subjeto
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % ERD = m_c/r_c
    ERD{y,1}= bsxfun(@times,m_c{y,1},1./r_c{y,1}) - 1;
end

%%
save(name,'ERD')
% ard = zeros();

%% Headplot for channels (lv = 1:n)

sub = subjeto;

hh = [1 3 5 7 9];
% hh = 1:12;
f = 1:length(hh);
da = zeros(length(hh),22,1750);
% ERD{sub}(11,7,:) = ERD{sub}(11,7,:)*-1;
for i =1:length(hh)
        da(i,:,:) = ERD{sub}(hh(i),:,:);
end

%% Headplot for all_channels
l1 = -5;
l2 = 5;
% -2.2359e+03 2.2359e+03
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];
i = 1;
 
for ch = 1:22
    figure(subjeto)
    subplot(6,7,posi(i))
    imagesc(t',f,squeeze(da(:,ch,:)),[l1 l2])
    title(['Ch: ' num2str(ch)])
    i = i+1;
end

% ard=MultiWP(X{1}{1},'sym2',4);