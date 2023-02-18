%% Load in the 

Dir1 = '/Users/bolger/Documents/work/Projects/Project-MorphemeB/';
filenom = 'MorphemeB_covariate_data.xlsx';

% Read line 1 for covariate type.
% Read lines 2 to 12853 and columns A to Q

[~,~,alldata] = xlsread(strcat(Dir1,filenom),'A1:Q12853');
types = alldata(1,:);
mots = alldata(:,ismember(types,'target_old'));
mots = mots(2:end,1);
lextype = alldata(:,ismember(types,'lex'));
lextype = lextype(2:end,1);
itemtype = alldata(:,ismember(types,'itemtype'));
itemtype = itemtype(2:end,1); 
afftype = alldata(:,ismember(types,'aff'));
afftype = afftype(2:end,1); 
wordtypes = cellfun(@strcat,itemtype,lextype,'UniformOutput',false);
wordtypes = cellfun(@strcat,wordtypes,afftype,'UniformOutput',false); 
wordtypes_ind = unique(wordtypes); 

%% Divide the wordtypes into the 9 different categories.

WordCond = cell(length(wordtypes_ind),1);

for counter = 1:length(wordtypes_ind)
    
    i = ismember(wordtypes,wordtypes_ind{counter,1});
    x = mots(i);
    WordCond{counter,1} = unique(x);
    
end
clear counter; 
%% Load in the *.wav files for each experimental condition. 

dirwav = '/Volumes/deepassport/Projects/Project_MorphemB/Experiment/';
D = dir(strcat(dirwav,'*.wav'));
fnoms = {D.name};
Fnames = {fnoms{~[D.isdir]}};
Dur = cell(length(WordCond),1);
Duravg = zeros(length(WordCond),1); 
Durtotal = zeros(length(WordCond),1);
Durstd = zeros(length(WordCond),1);

for counter = 1:length(WordCond)
    i2 = []; 
    snddur = zeros(length(WordCond{counter,1}),1);
    for cnt = 1:length(WordCond{counter,1})
        i2 = cellfun(@isempty,(strfind(Fnames,WordCond{counter,1}(cnt))),'UniformOutput',false);
        I = [cell2mat(i2) == 0];
        [sndin, fs] = audioread(strcat(dirwav,char(Fnames(I))));
        snddur(cnt) = length(sndin)/fs;
    end 
    Dur{counter,1} = snddur.*1000; 
    Durtotal(counter) = sum(snddur)*1000;
    Duravg(counter) = mean(snddur)*1000; 
    Durstd(counter) = std(Duravg(counter));
end

Dur2 = reshape(cell2mat(Dur),[length(WordCond{counter,1}) length(WordCond)]);

figure;
boxplot(Dur2,'notch','on','Labels',wordtypes_ind)

figure;
boxplot(Dur2,'notch','on','Labels',{'suffixed nonword-NS' 'non-suffixed nonword-NS' 'word-NS' 'suffixed nonword-PS' ,...
    'non-suffixed nonword-PS' 'word-PS' 'suffixed nonword-TS' 'non-suffixed nonword-TS' 'word-TS' })

figure;
boxplot(cell2mat(Dur),'notch','on','Labels','All Conditions')

[fid, p] = fopen(strcat(Dir1,'item_meandur.txt'),'w+');
fprintf(fid,'%d\n',Duravg);
fclose(fid); 

