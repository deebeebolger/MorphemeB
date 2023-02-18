% EEGLAB history file generated on the 12-Sep-2018
% ------------------------------------------------
%% LOAD SET FILES INTO EELAB SESSION

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; 

currdirectory = fullfile('Users','bolger','Desktop','Files-T0Adjusted',filesep);
savedirectory = fullfile(currdirectory,'blc',filesep);

files_cond= dir(currdirectory);
currfiles=files_cond;
fileIndex= find(~[currfiles.isdir]);
FInd=fileIndex;
filenum=dir(strcat(currdirectory,'*.set'));
filenom={filenum.name};                           %titles of the all the .set files to be loaded
EEG = pop_loadset('filename',filenom,'filepath',currdirectory);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
EEG=eeg_checkset(EEG); eeglab redraw;


for counter = 1:length(ALLEEG)
    
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',counter,'study',0);
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
    EEG = pop_rmbase( EEG, [-200 0]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    
    name_new = strcat(ALLEEG(counter).setname,'-bl');
    EEG = pop_saveset( EEG, 'filename',name_new,'filepath',savedirectory);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    eeglab redraw;
    
end
