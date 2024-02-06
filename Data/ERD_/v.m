

data = 1:999; %Dummy data. Just for testing.

Fs = 8000; % All the songs we'll be working on will be sampled at an 8KHz rate

tWindow = 64e-3; % The window must be long enough to get 64ms of the signal
NWindow = Fs*tWindow; % Number of elements the window must have
window = hamming(NWindow); % Window used in the spectrogram

NFFT = 512;
NOverlap = NWindow/2; % We want a 50% overlap

[S, F, T] = spectrogram(data, window, NOverlap, NFFT, Fs);