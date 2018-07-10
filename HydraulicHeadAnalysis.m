%% For the analysis load the file HEADS.mat
% This is not included in the repository because the size is about 224 MB
%
clear TOPHEAD
% For each time step compute the top most head value
for it = 1:size(HEADS,1)
    it
    TOPHEAD{it,1} = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
    for ii = 1:size(HEADS{1,1},1)
        for jj = 1:size(HEADS{1,1},2)
            for kk = 1:size(HEADS{1,1},3)%:-1:1
                if HEADS{it,1}(ii,jj,kk) ~= HEADS{1,1}(1,1,1)
                    TOPHEAD{it,1}(ii,jj) = HEADS{it,1}(ii,jj,kk);
                end
            end
        end
    end
end
%
%% Calculate Standard deviation
head_std = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
for ii = 1:size(HEADS{1,1},1)
    for jj = 1:size(HEADS{1,1},2)
        temp = nan(size(TOPHEAD,1),1);
        for kk = 1:1:size(TOPHEAD,1)
            temp(kk,1) = TOPHEAD{kk,1}(ii,jj);
        end
        head_std(ii,jj) = std(temp);
    end
end
%%
temp = head_std;
temp(head_std>2) = nan;
surf(temp,'edgecolor','none')
view(360,-90)