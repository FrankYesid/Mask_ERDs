function EEG_ = struct_EEGLAB(data,channels);
load('EEG.mat')
EEG_ = EEG;
EEG_.nbchan        = numel(channels);
EEG_.trials            = size(data,1);
EEG_.pnts             = size(data,3);
EEG_.srate             = 1000;
EEG_.xmin             = -200;
EEG_.xmax            = 800;
EEG.times              = -200:800-1/1000;
EEG_.data             = permute(data(:,channels,:),[2,3,1]);

EEG_.icaact          = [];
EEG_.icawinv        = [];
EEG_.icasphere    = [];
EEG_.icaweights   = [];
EEG_.icaweights   = [];
EEG_.icachansind = [];

EEG_.nbchan = numel(channels);
EEG_.data = EEG_.data(channels,:,:);

end
