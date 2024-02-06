clear
load('F:\Database\BCICIV_2a\BCICIV_2a.mat')
class = 1;
for s = [3,8,9]
    d = X{s}(y{s}==class);
    for tr = 1:numel(d)
        dats{tr} = d{tr}.^2;
    end
    for fil = 1:35
        dat{s}{fil} = fcnfiltband(dats,fs,[fil+3 fil+5],5); % suj x fil x trials x (time x ch)
    end
end

for s = [3,8,9]
    for fil = 1:35
        for tr = 1:numel(dat{s}{fil})
            for ch = 1:22
                X_real{s}{fil}(ch,tr,:) = dat{s}{fil}{tr}(:,ch);        % suj x ch x (fil x trial x time)
            end
        end
    end
end


%% con teager

% load('F:\Database\BCICIV_2a\BCICIV_2a.mat')
class = 1;
order = 3;
framelen = 11;

for s = [3]
    d = X{s}(y{s}==class);
    for tr = 1:numel(d)
        aa{tr} = sgolayfilt(d{tr},order,framelen);
    end
    clear d
    d = aa;
    for tr = 1:numel(d)
        dats1{tr} = d{tr}.^2;
    end
    for tr = 1:numel(d)
        for ch = 1:22
            da(:,ch) = abs(teager_2(d{tr}(:,ch)));
        end
        dats2{tr} = da.^2;
    end
    for tr = 1:numel(d)
        for ch = 1:22
%             d2 = mean(dats1{tr}(1748:1750,ch));
%             dats1{tr}(1748,ch) = d2;
            xpru{s}(tr,:,ch) = bsxfun(@times,dats1{tr}(1:1748,ch),1./dats2{tr}(:,ch)) - 1;
        end        
    end
%     for fil = 1:35
%         dat{s}{fil} = fcnfiltband(dats,fs,[fil+3 fil+5],5); % suj x fil x trials x (time x ch)
%     end
end

%%
for s = [3]
%     for fil = 1:35
        for tr = 1:numel(X{s})%(dat{s}{fil})
            for ch = 1:22
%                 X_suj{s}{fil}(ch,tr,:) = dat{s}{fil}{tr}(:,ch);        % suj x ch x (fil x trial x time)
                X_{s}{ch}(:,tr) = X{s}{tr}(:,ch);        % suj x ch x (fil x trial x time)
            end
        end
%     end
end

%% ERD-normal
ta = 1;
tb = 1*fs;
for  s = [3]
    for ch =1:35
        r_c = squeeze(mean(X_real{s}{ch}(:,:,ta:tb),3));
        r_ = squeeze(mean(r_c,2));
        m_c = squeeze(mean(X_real{s}{ch},2));
        ERD{s}{ch} =  (bsxfun(@times,m_c,1./r_) - 1);
    end
end

%% ERD-normal con teager
ta = 1;
tb = 1*fs;
for  s = [3]
    for ch =1:35
        r_c = squeeze(mean(X_suj{s}{ch}(:,:,ta:tb),3));
        r_ = squeeze(mean(r_c,2));
        m_c = squeeze(mean(X_real{s}{ch},2));
        ERD{s}{ch} =  (bsxfun(@times,m_c,1./r_) - 1);
    end
end

%% plot errorbar
set(0,'DefaultFigureWindowStyle','docked')
for s = [3]
    figure
    for ch = 1:22
        hold on
%         errorbar(mean(ERD{s}{ch}),std(ERD{s}{ch}))
        plot(mean(ERD{s}{ch}))
        xticks([250 500 750 1000 1250 1500 1750])
        xticklabels({'1','2','3','4','5','6','7'})
        xlim([1 1750])
    end
    title(['Sujeto ' num2str(s)])
end

%%
