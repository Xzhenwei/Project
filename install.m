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

epsilon=1e-2;
SeismicEulerWorkBook(epsilon)
char=['Epsilon',num2str(epsilon),'.mat'];
save(char,'-mat')

epsilon=20;
SeismicEulerWorkBook(epsilon)
char=['Epsilon',num2str(epsilon),'.mat'];
save(char,'-mat')

epsilon=30;
SeismicEulerWorkBook(epsilon)
char=['Epsilon',num2str(epsilon),'.mat'];
save(char,'-mat')
epsilon=45;
SeismicEulerWorkBook(epsilon)
char=['Epsilon',num2str(epsilon),'.mat'];
save(char,'-mat')
end


