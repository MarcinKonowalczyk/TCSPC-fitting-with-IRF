close all; clear; clc;

%% Find all the .m files
all_files = dir(pwd);
m_files = {};
for j = 1:numel(all_files)
    f = all_files(j).name;
    if numel(f)>2 && strcmp(f(end-1:end),'.m')
        if ~strcmp(f(1:end-2),mfilename)
            m_files(end+1) = {f}; %#ok<SAGROW>
        end
    end
end

%% Check the code dependencies
for j = 1:numel(m_files)
    filename = m_files{j};
    fprintf('%s\n',filename);
    [flist,plist] = matlab.codetools.requiredFilesAndProducts(filename);
    fprintf('Dependencies:\n');
    for k = 1:numel(flist)
        fprintf(' %s\n',(flist{k}));
    end
    fprintf('Toolboxes:\n');
    for k = 1:numel(plist)
        fprintf(' %s\n',(plist(k).Name));
    end
    fprintf('\n');
end
