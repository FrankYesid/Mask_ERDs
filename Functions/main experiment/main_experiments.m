%% main_experiments.m
%   El presente script permite clasificar una base de datos de BCI respecto
%   a diferentes metodología propuestas.
% 

clear all; close all;clc;warning off

% Cargar base de datos
% load BCICIV_1.mat
load BCICIV_2a.mat
% load objcv


% Define working labels
labels = [1 2]; % 

% Elegir tipo de experimento
experiment = 2;

% Dependiendo de el experimento se ejecuta el c�digo respectivo.
switch experiment
    case 1 %CSP
        fprintf('Running CSP\n')
        CS_P
    case 2
        fprintf('Running ERDSCSP\n')
        ERDSCSP
    otherwise
end
