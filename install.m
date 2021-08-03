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

% addpath(fullfile(maindir, 'QuarterCar'));
addpath(fullfile(maindir, 'Seismic'));
% addpath(fullfile(maindir, 'vonKarmanBeam'));

run SeismicWorkBook2.mlx
wsTime = clock;
char=['Date_',num2str(wsTime(2)),num2str(wsTime(3)),...
    ' Time_',num2str(wsTime(4)),num2str(wsTime(5)),'.mat'];
save(char,'-mat')
end


