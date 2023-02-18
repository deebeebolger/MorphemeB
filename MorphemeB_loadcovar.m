
Dir1 = '/Users/bolger/Documents/work/Projects/Project-MorphemeB/';
filenom = 'MorphemeB_covariate_data.xlsx';

% Read line 1 for covariate type.
% Read lines 2 to 12853 and columns A to Q

[~,~,alldata] = xlsread(strcat(Dir1,filenom),'A1:Q12853');
types = alldata(1,:);
mots = alldata(:,ismember(types,'target_old'));
mots = mots(2:end,1);
upoint = alldata(:,ismember(types,'UPphon'));
upoint = upoint(2:end,1);
wfreq = alldata(:,ismember(types,'WFreq'));
wfreq = wfreq(2:end,1);
spfreq = alldata(:,ismember(types,'SPFreq'));
spfreq = spfreq(2:end,1);
nletter = alldata(:,ismember(types,'Letters'));
nletter = nletter(2:end,1);
lextype = alldata(:,ismember(types,'lex'));
lextype = lextype(2:end,1);

%% Load in the subject-level epoched data.



[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;               
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

eegdir = '/Users/bolger/Documents/work/Projects/Project-MorphemeB/EEG_Data/';
wordtype = 'NonWord_S';
wordcond = 'NS';
fulldir = strcat(eegdir,wordtype,'/',wordcond,'/');
D = dir(strcat(fulldir,'*.set'));
fnoms = {D.name};
Fnames = {fnoms{~[D.isdir]}};

EEG = pop_loadset('filename',Fnames,'filepath',fulldir);
EEG = eeg_checkset( EEG );
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
eeglab redraw; 

%% %  For each word (mot) find instances in the events field and add the
%  uniqueness point (upoint), the wfreq, the spfreq and the number of
%  letter. Record these in the events file. 


for counter = 1:size(ALLEEG,2)
    
    Ecurr = ALLEEG(1,counter).event;
    motscurr = {Ecurr.type};
    upphon = zeros(length(motscurr),1);
    wordfreq = zeros(length(motscurr),1);
    spkfreq = zeros(length(motscurr),1);
    letters = zeros(length(motscurr),1);
    
    for i = 1:length(motscurr)
        
        indx = cellfun(@isempty,strfind(mots,motscurr{1,i}(3:end)),'UniformOutput',false);
        indx = cell2mat(indx);
        Indx = find(indx==0);
        Indx = Indx(1);
        letters(i,1) = nletter{Indx,1};
        
        if strcmp(lextype{Indx,1},'word') && ~strcmp(upoint{Indx,1},'#N/A')
            upphon(i,1) = upoint{Indx,1};
            wordfreq(i,1) = str2double(wfreq{Indx,1});
            spkfreq(i,1) = str2double(spfreq{Indx,1});
        else
            disp('NonWord');
        end
     ALLEEG(1,counter).event(i).UPphon = upphon(i,1);
     ALLEEG(1,counter).event(i).WFreq = wordfreq(i,1);
     ALLEEG(1,counter).event(i).SPFreq = spkfreq(i,1);
     ALLEEG(1,counter).event(i).Letters = letters(i,1);     
    end
    
     EEG = pop_saveset(ALLEEG(1,counter), 'filename',ALLEEG(1,counter).setname,'filepath',fulldir);  % Saves a copy of the current resampled dataset to the current directory
     eeglab redraw
    
end









