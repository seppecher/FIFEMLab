s = what('fifemlab');
disp('fifemlab setup');
addpath(s.path);
addpath([s.path '/classes']);
addpath([s.path '/data']);
addpath([s.path '/documentation']);
addpath([s.path '/functions']);
savepath;
disp('setup done');