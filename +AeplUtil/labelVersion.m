function labelVersion(path)

% GIT = tdfread('.git/logs/HEAD');
% disp(GIT)
% Y = GIT.clone0x3A_from_https0x3A0x2F0x2Fgithub0x2Ecom0x2Fdeep0x2Dfried0;
% X = GIT.x0000000000000000000000000000000000000000_542fa38aff1b597ceb0fa;
% 
% x = X(end,:);
% y = Y(end,:);

[~,git_hash_string] = system('git rev-parse HEAD');
[~,git_log] = system('git --no-pager log --oneline -1');

fileID = fopen(fullfile(path,'Commit Stamp.txt'),'w');
fprintf(fileID,git_hash_string);
fprintf(fileID,'    ');
fprintf(fileID,git_log);

fclose(fileID);
end 

