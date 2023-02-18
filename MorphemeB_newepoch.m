%% FIND INCORRECT AND LATE RESPONSES
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;                %open eeglab session
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

rincorr = '252';        %define trigger for incorrect response
rcorr = '251';          %define the trigger for correct response
delay = '100';

sujnum = [1:27];
dirbase = '/Users/bolger/Documents/work/Projects/Project-MorphemeB/EEG_Data/';

for scount = 1:length(sujnum)
    
    %% Construct current subject-level directory.
    
    sujcurr = [];
    if scount<10 && scount>1
        sujcurr = strcat('S0',num2str(scount));
    elseif scount>=10
        sujcurr = strcat('S',num2str(scount));
    elseif scount ==1
        sujcurr = 'P01';
    end
    
    dirsuj = strcat(dirbase,sujcurr,'/');
    
    % Only load in the merged files.
    file2load = dir(strcat(dirsuj,'*merged.set'));
    EEG = pop_loadset('filename',file2load.name,'filepath',dirsuj);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
    eeglab redraw;
    
    %% Find incorrect and slow answers and mark. 
    
    Incorr_indx = find(ismember({EEG.event.type},rincorr));
    Corr_indx = find(ismember({EEG.event.type},rcorr));
    delays = find(ismember({EEG.event.type}, delay));
    
    if isempty(Incorr_indx)
        
        display('Excellent! Aucune faute!')
        
    elseif ~isempty(Incorr_indx)
        
        display( strcat(num2str(length(Incorr_indx)),' fautes! Il faut les trouver!'))
        for evnt_counter=1:length(Incorr_indx)
            
            EEG.event(Incorr_indx(evnt_counter)-1).type = '999';
        end
        
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = eeg_checkset( EEG );
        
    end
    
    
    
    if isempty(delays)
        
        display('Excellent! No slow responses!')
        
    elseif ~isempty(delays)
        
        display( strcat(num2str(length(delays)),' fautes! Il faut les trouver!'))
        for evnt_counter = 1:length(delays)
            
            EEG.event(delays(evnt_counter)-2).type = '998';
        end
        
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = eeg_checkset( EEG );
        
    end
    
    
    EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
    eeglab redraw
    
    beforeEp=CURRENTSET;   % define the index of the dataset to be segmented
    
    %% Segment Data Continuous Data.
    % Begin by segmenting all trials based on word data only.
    trial_lims = [-0.2 1];
    unwanted_trigs = {'0' '100' '251' '252' 'boundary' '999' '998'};
    condunique = unique({EEG.event.type});
    unwanted_i = logical(~ismember(condunique,unwanted_trigs));
    condsoi = condunique(unwanted_i);
    
    fnom = strcat(EEG.setname,'-epochs-allwords');
    EEG = pop_epoch( EEG, condsoi, [trial_lims(1) trial_lims(2)], 'newname', char(fnom), 'epochinfo', 'yes');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom),'gui','off');
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
    eeglab redraw
    
    %% Baseline Correct the segmented data.
    
    fnombl = strcat(fnom,'-bl');
    EEG = pop_rmbase( EEG, [trial_lims(1)*1000 0]);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnombl),'filepath',EEG.filepath);
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
end




