% Add SOSTOOLS to path
%

cm = computer;
if cm(1) == 'M' ||  cm(1)=='G'
    addpath(pwd);
    addpath(genpath([pwd '/multipoly']));
    addpath(genpath([pwd '/internal']));
    addpath(genpath([pwd '/demos']));
    addpath(genpath([pwd '/custom']));
    addpath(genpath([pwd '/dpvar']));
elseif cm(1) == 'P'
    addpath(pwd);
    addpath(genpath([pwd '\multipoly']));
    addpath(genpath([pwd '\internal']));
    addpath(genpath([pwd '\demos']));
    addpath(genpath([pwd '\custom']));
    addpath(genpath([pwd '\dpvar']));
end

% Alternative single-line syntax to add all subfolders
% addpath(genpath(pwd));    