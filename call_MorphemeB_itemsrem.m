clear all;
close all;

path2items = '/Users/bolger/Documents/work/Projects/Project-MorphemeB/words-to-reject.txt';
Dirbase = fullfile('Users','bolger','Desktop','Files-T0Adjusted',filesep);
sujnoms = {'s16' 's17' 's18' 's20' 's21' 's22' 's23' 's24' 's25' 's26' 's27' 's28'};

for count = 1:length(sujnoms)
    items_indx = MorphemeB_itemsrem(sujnoms{count},Dirbase,path2items);
end