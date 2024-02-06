%% Database
load('BCICIV_2a.mat')
%% Save Database of 1 subject
for i = 1:9
    x = X{i};
    tam = size(x,1);
    da = zeros(size(x{1},2),size(x{1},1),tam);
    for j = 1:tam
        da(:,:,j) = x{j}';
    end
    save(['BCICIV_2a' num2str(i) '.mat'],'da')
end
%% Database of 1 class
% Database
load('BCICIV_2a.mat')
clc;
class = 1;
for i = 1:9
    x = X{i}(y{i}==class);
    tam = size(x,1);
    da = zeros(size(x{1},2),size(x{1},1),tam);
    fprintf(['sujeto: ' num2str(i) ' de 9 \n'])
    for j = 1:tam
        da(:,:,j) = x{j}';
    end
    n{i,1} = da;
end
% ERD para toolbox
tiem = 250;
t1 = 1;
t2 = tiem*2;
for i = 1:9
    rnc{i,1} = squeeze(mean(n{i}(:,t1:t2,:),2));
end
for i = 1:9
    mc{i,1} = permute(n{i},[1 3 2]);
end
for i = 1:9
    ERD{i,1}= bsxfun(@times,mc{i},1./rnc{i,1}) - 1;     % channel,trial,time.
end
for i = 1:9
    ERDs{i,1} = permute(ERD{i},[1 3 2]);                % channel,time,trial.
end
%% guardar erds para cada sujeto 
for i = 1:9
     tam = ERDs{i};
    fprintf(['sujeto: ' num2str(i) ' de 9 \n'])
    save(['BCICIV_2a' num2str(i) 'ERD'],'tam')
end