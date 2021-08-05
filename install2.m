function install2
% run this file first to install all external packages  and
% also the main software appropriately

maindir = fileparts(mfilename('fullpath'));
extdir = 'ext';
cocoinstall = fullfile(maindir, extdir , 'coco','startup.m');
run(cocoinstall);

addpath(fullfile(maindir, extdir,'combinator'));

addpath(fullfile(maindir, extdir, 'tensor_toolbox'));

addpath(genpath(fullfile(maindir, extdir, 'YetAnotherFEcode','src')));

addpath(genpath(fullfile(maindir, extdir, 'YetAnotherFEcode','examples')));

addpath(fullfile(maindir, extdir, 'Wrappers'));

addpath(fullfile(maindir, 'src'));

addpath(fullfile(maindir, 'src','misc'));

addpath(fullfile(maindir, 'src','multiindex')); 

addpath(fullfile(maindir, 'src', 'frc'));

% addpath(fullfile(maindir, 'QuarterCar'));
addpath(fullfile(maindir, 'vonKarmanBeam'));

epsilon=1e-2;
vonKarmanEuler(epsilon)
clear all
epsilon=1.5e-2;
vonKarmanEuler(epsilon)
clear all
epsilon=2e-2;
vonKarmanEuler(epsilon)
clear all
epsilon=3e-2;
vonKarmanEuler(epsilon)

end


