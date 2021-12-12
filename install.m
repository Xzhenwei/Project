function install
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

addpath(fullfile(maindir, 'src', 'psd'));
cd('QuarterCar')
run QuarterCarWorkbook.mlx
cd('..')
cd('Seismic')
run SeismicWorkBook2.mlx  
cd('..')
cd('vonKarmanBeam')
run vonKarmanBeamWorkBook.mlx 
cd('..')
cd('vonKarmanShell')
run vonKarmanShellWorkbook.mlx
cd('..')

% cd('Seismic')
% cd('..')
% cd('vonKarmanBeam')
% addpath(fullfile(maindir, 'QuarterCar'));
% addpath(fullfile(maindir, 'Seismic'));
% addpath(fullfile(maindir, 'vonKarmanBeam'));
% addpath(fullfile(maindir, 'vonKarmanShell'));
% psd_full (0.8)
% run SeismicWorkBook2.mlx
% save('linearSeismic.mat')
end


