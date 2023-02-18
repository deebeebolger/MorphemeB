dirbase = fullfile('Users','bolger','Desktop','Files-T0Adjusted',filesep);

Group = {'Word' 'NonWord-S' 'NonWord-NS'};
Groupcode = [10 20 30];
SubGroup = {'TS' 'PS' 'NS'};
SubGroupcode = [1 2 3];

for scount = 1:length(ALLEEG)
    
   [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',scount,'study',0);
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
    for counter = 1:length(Group)
        for counter2 = 1:length(SubGroup)
            
            Gcond = Group{1,counter};
            SGcond = SubGroup{1,counter2};
            Gcondcode = Groupcode(strcmp(Gcond,Group));
            SGcondcode = SubGroupcode(strcmp(SGcond,SubGroup));
            cond = num2str(Gcondcode + SGcondcode);
            
            disp(horzcat('Current condition is ',Gcond,'-',SGcond,': Condition code ',cond));
            
            %% Need to find instances of this condition in the current events field of the EEG structure.
            
            currevents = {EEG.event.type};
            X = strfind(currevents,cond);
            x= cell2mat(cellfun(@isempty,X,'UniformOutput',false));
            CondCurr = {EEG.event(x==0).type};   % The condition names that will be used for segmentation.
            
            
            
            fnom = strcat(EEG.setname(1:4),Gcond,'-',SGcond);
            savepath = fullfile(dirbase,[Gcond,'-',SGcond],filesep);
            EEG = pop_selectevent( EEG, 'type',CondCurr,'deleteevents','on','deleteepochs','on','invertepochs','off');
            EEG.setname= fnom;
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom),'overwrite','off'); % current set = xx;
            EEG = eeg_checkset( EEG );
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            eeglab redraw
            EEG = pop_saveset( EEG, 'filename',char(fnom),'filepath',savepath);  % Saves a copy of the current resampled dataset to the current directory
            eeglab redraw
            
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',scount,'study',0);
            EEG = eeg_checkset( EEG );
            eeglab redraw
            
        end
    end
end







