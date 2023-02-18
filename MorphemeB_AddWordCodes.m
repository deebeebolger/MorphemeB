%% Function to read the individual word data from excel files and integrate into the current EEG data structure.
% The excel file sheet that should be used differs depending on the subject
% number.
% For each subject, the three excel files corresponding to the three lists
% should be called (RandomiseList1.xlsx, RandomiseList2.xlsx and
% RandomiseList3.xlsx). 
%***********************************************************************************************************

% Load in the trial-level information from the excel files.
% Define the path to the excel files. 
Direxcel = '/Users/bolger/Documents/work/Projects/Project-MorphemeB/';
xlfile_List = {'RandomiseList1.xlsx' 'RandomiseList2.xlsx' 'RandomiseList3.xlsx'};
sheetnames = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
ListData = [];          % Define an empty structure.
ListData.List1.dir = strcat(Direxcel,xlfile_List{1,1});
ListData.List2.dir = strcat(Direxcel,xlfile_List{1,2});
ListData.List3.dir = strcat(Direxcel,xlfile_List{1,3});

itemtype = cell(3,length(sheetnames)); 
lex = cell(3,length(sheetnames));
affix = cell(3,length(sheetnames));
codes = cell(3,length(sheetnames));
words = cell(3,length(sheetnames));

for lcount = 1:3
    currlist = strcat('List',num2str(lcount));
    for sheetnum = 1:length(sheetnames)
        [~,itemtype{lcount,sheetnum},~] = xlsread(ListData.(matlab.lang.makeValidName(currlist)).dir,sheetnames{1,sheetnum},'E2:E154');
        [~,lex{lcount,sheetnum},~] = xlsread(ListData.(matlab.lang.makeValidName(currlist)).dir, sheetnames{1,sheetnum},'G2:G154');
        [~,affix{lcount,sheetnum},~] = xlsread(ListData.(matlab.lang.makeValidName(currlist)).dir, sheetnames{1,sheetnum},'H2:H154');
        [codes{lcount,sheetnum},~,~] = xlsread(ListData.(matlab.lang.makeValidName(currlist)).dir, sheetnames{1,sheetnum},'K2:K154');
        [~,words{lcount,sheetnum},~] = xlsread(ListData.(matlab.lang.makeValidName(currlist)).dir, sheetnames{1,sheetnum},'L2:L154');
        
        
    end
    ListData.(matlab.lang.makeValidName(currlist)).itemtype = [itemtype{lcount,:}];
    ListData.(matlab.lang.makeValidName(currlist)).lexicality = [lex{lcount,:}];
    ListData.(matlab.lang.makeValidName(currlist)).affixed = [affix{lcount,:}];
    ListData.(matlab.lang.makeValidName(currlist)).codes = [codes{lcount,:}];
    ListData.(matlab.lang.makeValidName(currlist)).words = [words{lcount,:}];
end

sheetsubj = [1 2 15 22;3 4 16 23;5 6 17 24; 7 8 18 25;9 10 19 26; 11 12 20 27; 13 14 21 28];   % Read in rows.
Dirbase = '/Users/bolger/Documents/work/Projects/Project-MorphemeB/EEG_Data/';      % The default path.

[ALLEEG EEG CURRENTSET, ALLCOM] = eeglab;    %open eeglab session

Listcurr = 'List3';
for scount = 1:size(sheetsubj,1)                    % For each sheet
    for fcount = 2:size(sheetsubj,2)                % For each subject group
        
        dircurr= [];  
        if sheetsubj(scount,fcount)==1% Define the current folder.
            dircurr = strcat(Dirbase,'P0',num2str(sheetsubj(scount,fcount)),'/');
        elseif sheetsubj(scount,fcount)<10 && sheetsubj(scount,fcount)>1
            dircurr = strcat(Dirbase,'S0',num2str(sheetsubj(scount,fcount)),'/');
        elseif sheetsubj(scount,fcount)>=10
            dircurr = strcat(Dirbase,'S',num2str(sheetsubj(scount,fcount)),'/');
        end
        
        %% Load in the current subject. 
        D = dir(strcat(dircurr,'*.set'));
        fnoms = {D.name};
        Fnames = {fnoms{~[D.isdir]}};
        % Find the names corresponding to the current list (List1, List2 or
        % List3).
        indxl = cell2mat(cellfun(@isempty,strfind(Fnames,Listcurr),'UniformOutput',false));
        Fnames_curr = {Fnames{[indxl==0]}}; 
        
        
        EEG = pop_loadset('filename',Fnames_curr,'filepath',dircurr);
        EEG = eeg_checkset( EEG );
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
        eeglab redraw; 
        
        %% Convert from string to number if in string format
        if ischar(EEG.event(2).type)
            for ecnt1 = 1:length(EEG.event)
                EEG.event(ecnt1).type = str2double(EEG.event(ecnt1).type);
            end
        else
            disp('Not in string format, no need to change!');
        end
        %%  Call of function to correct the triggers.
        trigcorr = Trigger_Fixit([EEG.event.type]);  
        for ecnt = 1:length(trigcorr)
            EEG.event(ecnt).type = trigcorr(ecnt);
        end
        EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',dircurr);  % Saves a copy of the current resampled dataset to the current directory
        eeglab redraw
        
        %% 
        words_curr = [];
        words_curr = {ListData.(matlab.lang.makeValidName(Listcurr)).words{:,scount}};
        codeu = unique(ListData.(matlab.lang.makeValidName(Listcurr)).codes(:,scount));
        wcnt = 1;
        
        for ecnt = 1:length(trigcorr)
            x = EEG.event(ecnt).type;
            if ismember(x,codeu)
            wcnt
            ibidon = strfind(words_curr{1,wcnt},'_');    
            EEG.event(ecnt).type = strcat(num2str(x),words_curr{1,wcnt}(ibidon+1:end-4));
            wcnt = wcnt +1; 
            else
                EEG.event(ecnt).type = num2str(x);
            end
            
        end
        EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',dircurr);  % Saves a copy of the current resampled dataset to the current directory
        eeglab redraw
        
        
        
        
    end
end


