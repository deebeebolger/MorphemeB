% Date: 01-02-1976   Programmed by: Deirdre Bolger
% Script to carry out a non-parametric statistic test on ERP data. It uses
% a cluster-based test statistic to deal with the multiple comparisons
% problem (MCP).
% The main steps of the test are as follows :
% 1. The experimental conditions are compared for every channel-time pair by
% means of a t-test.
% 2. All those samples whos t-value exceeds a predefined
% threshold are selected.
% 3.The selected samples are clustered in connected sets on the basis of temporal and spatial adjacency.
% 4. The sume of the t-values within each cluster are taken which yields
% the cluster-level statistics.
% 5. The maximum of the cluster-level statistics are found, which gives the
% test statistic used to evaluate the effect of the experimental
% conditions.
% As the MCP is resolved by the cluster-based test statistic, the
% significance probability has to be calculated using the Monte Carlo
% method.
% Very important: Need to know if the structure of your experiment falls
% under the heading of a "between-" or a "within" unit of observation
% design. The Monte Carlo significance probability is calculated
% differently depending on the experimental structure (within or between UO
% design).

clear all;
close all;


%% Define the condition names and groups and root directory *************************
Groups= {};
Conds = {};
condnoms1=strcat(Groups{1,1},'-',Conds{1,1});
condnoms2=strcat(Groups{1,2},'-',Conds{1,2});
condnames=cell(1,2);
condnames{1,1}=condnoms1;
condnames{1,2}=condnoms2;
dir_root = fullfile(filesep,'',filesep); % %%BESOIN DE CHANGER CE CHEMIN!!!!
Time_Interval = []; %time interval over which to carry out the analysis in seconds
timestep=0.05;
SLOption='No';         %choose to calculate the surface laplacian
channum=64;
layout_dir = 'tobeadded/MATLAB/fieldtrip-master/template/layout/biosemi64.lay' ; %% BESOIN DE CHANGE CE CHEMIN!!!!
anal_type = 'multielec'; %single_elec ou multielec
timestep_anal ='yes'; %or "yes"
filterdata ='no'; %or 'non'

%% **************************************Load in the data files and convert to Fieldtrip form for spatial permutation analysis*****************************************

if strcmp(Groups{1,1},Groups{1,2})==1
    lenGroup=1;
else
    lenGroup=length(Groups);
end

if strcmp(Conds{1,1},Conds{1,2})==1
    lenCond=1;
else
    lenCond=length(Conds);
end


AllFiles=cell(lenGroup,lenCond);
AllDirs=cell(lenGroup,lenCond);
DataIn_eeglab=cell(lenGroup,lenCond);
%DataIn_eeglab=cell(lenGroup,size(Conds,2));

for Gcounter=1:lenGroup
    for Condcounter=1:lenCond
        
        [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
        
        currfolder = fullfile(Groups{1,Gcounter},filesep,Conds{1,Condcounter});
        AllDirs{Gcounter,Condcounter} = fullfile(dir_root,currfolder,filesep);
        AllFiles{Gcounter,Condcounter}=dir(strcat(AllDirs{Gcounter,Condcounter},'*.set'));   %declare empty variable for files
        currfiles=AllFiles{Gcounter,Condcounter};
        filenom={AllFiles{Gcounter,Condcounter}.name};
        
        EEG=pop_loadset('filename',filenom,'filepath',AllDirs{Gcounter,Condcounter});
        [ALLEEG, EEG, CURRENTSET]=pop_newset(ALLEEG,EEG,0,'study',0);
        EEG=eeg_checkset(EEG); eeglab redraw;
        
        DataIn_eeglab{Gcounter,Condcounter}=ALLEEG;
        numsuj=length(ALLEEG);
        
        
        if strcmp(filterdata,'yes')
            display('*********************************Lowpass filtering using a FIR windowed sinc filter***********************************')
            SR=ALLEEG(1).srate;
            [M, dev]=pop_firwsord('blackman',SR, 2);    %calculate the filter order
            for eegcounter=1:numsuj
                display(ALLEEG(eegcounter).setname)
                [DataIn_eeglab{Gcounter,Condcounter}(eegcounter),~,~]=pop_firws(ALLEEG(eegcounter),'fcutoff',30,'forder',M,'ftype','lowpass','wtype','blackman');
            end
        end
        
    end
end

%% Check that both datasets (conditions) have the same T0.

if size(DataIn_eeglab,2)==1
    DataIn_eeglab=DataIn_eeglab';
end

[DataIneeg_new,Time]=CREx_alignT0(DataIn_eeglab);  % call of function CREx_alignT0 to realign datasets whose T0 differ
%DataIneeg_new=DataIneeg_new;

%% Convert from eeglab to fieldtrip

if lenGroup ==1
    DataIn=cell(numsuj,lenGroup,lenCond);
    DataIn_pp=cell(numsuj,lenGroup,lenCond);
elseif lenGroup >1
    DataIn = cell(numsuj,lenCond,lenGroup);
    DataIn_pp = cell(numsuj,lenCond,lenGroup);
    bidg = lenGroup;
    lenGroup = lenCond;
    lenCond = bidg;
end

for gcnt=1:lenGroup
    for ccnt=1:lenCond
        for fcnt=1:numsuj
            
            DataIn_pp{fcnt,gcnt,ccnt}= eeglab2fieldtrip(DataIneeg_new{gcnt,ccnt}(fcnt), 'preprocessing', []);
            
            alleegs=DataIneeg_new{gcnt,ccnt}(fcnt);
            cfg=[];
            cfg.channel={alleegs.chanlocs(1:64).labels};
            cfg.trials='all';
            cfg.keeptrials='no';
            cfg.covariance='no';
            
            DataIn{fcnt,gcnt, ccnt}=ft_timelockanalysis(cfg, DataIn_pp{fcnt,gcnt,ccnt});
        end
    end
end

if size(DataIn,2)==1
    DataIn=squeeze(DataIn);
end

%% OPTION TO CALCULATE THE SURFACE LAPLACIAN

if strcmp(SLOption,'Yes')==1
    
    disp('Calculating the Surface Laplacian')
    
    csdFileName='C:\Users\bolger\Documents\MATLAB_tools\CSDtoolbox\CSDtoolbox\resource\10-5-System_Mastoids_EGI129.csd';
    Labels=textread('C:\Users\bolger\Documents\MATLAB_tools\CSDtoolbox\CSDtoolbox\resource\E64_10_20.txt','%s');
    
    montage=CREx_ExtractMontage(csdFileName,Labels);    %Extract the Montage
    [Gout,Hout]=CREx_CalcGH(montage,7);                          %Calculate the G and H transformation matrices
    
    SL=cell(size(DataIn));
    CSD=cell(size(DataIn));
    
    for GCount=1:size(DataIn,2)
        for SCount=1:size(DataIn,1)
            
            [CSD{SCount,GCount},SL{SCount,GCount}]=CREx_CalcCSD(DataIn{SCount,GCount}.avg(1:64,:),Gout,Hout,64);
            DataIn{SCount,GCount}.avg=SL{SCount,GCount};
            %DataIn{SCount,GCount}.SL='yes';
            
        end
    end
end

%% Compute non-parametric statistical test


if lenGroup>1 && size(DataIn,3)>1
    DEV={DataIn{:,1,:}};
    STD={DataIn{:,2,:}};
    DEV_char=Conds{1,1};
    STD_char=Conds{1,2};
elseif lenGroup==1 && size(DataIn,3)==1
    DEV={DataIn{:,1}};
    STD={DataIn{:,2}};
    DEV_char=condnames{1,1};
    STD_char=condnames{1,2};
elseif lenGroup>1 && size(DataIn,3)==1
    DEV={DataIn{:,1}};
    STD={DataIn{:,2}};
    DEV_char=condnames{1,1};
    STD_char=condnames{1,2};
end



%%

if strcmp(anal_type,'multielec')
    %% CALCULATE THE CLUSTER-BASED NON-PARAMETRIC TEST
    
    % Calculate the neighbours for current electrode configuration
    cfg=[];
    cfg.layout = layout_dir;
    lcfg.rotate=90;
    layout1 = ft_prepare_layout(cfg,DataIn{1,1});
    figure
    ft_plot_lay(layout1,'label','yes','point','no','outline','yes','box','no');
    
    % Define the neighbours
    cfg=[];
    cfg.method= 'triangulation';
    %cfg.neighbourdist = 15;   %define the maximum distance between neighbouring sensors
    cfg.layout = layout1;
    cfg.channel={DataIn{1,1}.label{1:channum}};
    neighbs=ft_prepare_neighbours(cfg,DataIn{1,1});
    
    cfg=[];
    cfg.neighbours=neighbs;
    cfg.enableedit='no';
    cfg.layout=layout1;
    ft_neighbourplot(cfg)
    
    cfg=[];
    cfg.channel = {EEG.chanlocs(1:channum).labels};
    cfg.latency=Time_Interval;
    cfg.avgoverchan = 'no';
    cfg.avgovertime= 'no';
    cfg.parameter = 'avg';
    cfg.method = 'montecarlo';
    cfg.statistic = 'ft_statfun_depsamplesT'; % independent samples T-statistic to measure effect size at sample level.
    cfg.correctm = 'cluster';                          % method to deal with multiple comparisons
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.tail = 0;
    cfg.clustertail = 0;
    cfg.correcttail ='alpha';
    cfg.alpha = 0.025;
    cfg.minnbchan = 3;            %the minimum number of neighbourhood channels required for selected channel to be included in the clustering algorithm
    cfg.numrandomization = 2000;  %number of draws from the permutation distribution
    cfg.neighbours=neighbs;
    
    if lenGroup >1
        designLen = lenGroup;
    elseif lenGroup ==1
        designLen = lenCond;
    end
    
    design = zeros(2,size(DataIn,1)*size(DataIn,3)*designLen);  %initialise the design
    
    % Create the design matrix to store information about the UOs --for
    % MotInter...also applies for MorphemeB
    
    for gcount = 1:designLen
        starter=((gcount-1)*size(design,2)/2)+1;
        
        if gcount == 1
            design(1,starter:size(DataIn,1)) = 1:size(DataIn,1);
            design(2,starter:size(DataIn,1))= gcount;
        else
            design(1,starter:size(DataIn,1)*gcount) = 1:size(DataIn,1);
            design(2,starter:size(DataIn,1)*gcount) = gcount;
        end
        
    end
    
    clear starter starter1 gcount
    
    
    display(design);
    cfg.design=design;
    
    
    cfg.ivar=2;                   %specifies the row of the design matrix corresponding to the independent variable
    cfg.uvar=1;                  %specifies the row of the design matrix corresponding to the unit of observation ...here the subjects
    tstats_all=ft_timelockstatistics(cfg,DEV{:},STD{:});   %DEV{:},STD{:}
    Tstats = tstats_all;
    
    
    %% Calculate the difference between the two experimental conditions
    
    cfg=[];
    cfg.latency = [tstats_all.time(1) tstats_all.time(end)];
    avgDEV = ft_timelockgrandaverage(cfg,DataIn{:,1});   %, ;  DEV{:} %Calculate the grand average over all subjects
    avgSTD = ft_timelockgrandaverage(cfg,DataIn{:,2});  %, DataIn{:,4});STD{:}
    
    cfg=[];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    diff_DEVvsSTD_all=ft_math(cfg,avgDEV,avgSTD);
    
    %% Construct a boolean matrix indicating which sensors are members of the significant clusters
    % Construct two matrices: one for positive and one for negative clusters
    % Begin with the positive clusters
    ispos=isfield(tstats_all,'posclusters');
    if ispos==1
        if  isfield(tstats_all.posclusters,'prob')==1
            pos_cluster_tstatsall= [tstats_all.posclusters(:).prob];
            pos_sig=find(pos_cluster_tstatsall<= .025);
        else
            display('No significant positive clusters!')
            pos_sig=[];
        end
        
        if isempty(pos_sig)==0
            pos=ismember(tstats_all.posclusterslabelmat,pos_sig);      %construct a boolean matrix of which channel-time pairs are part of significant clusters
        elseif isempty(pos_sig)==1
            
            display('No significant positive clusters!')
        end
    elseif ispos==0
        
        display('No positive clusters!');
        pos_sig=[];
    end
    
    
    %Now do the negative clusters
    isneg=isfield(tstats_all,'negclusters');
    if isneg==1
        if  isfield(tstats_all.negclusters,'prob')==1
            neg_cluster_tstatsall=[tstats_all.negclusters(:).prob];
            neg_sig=find(neg_cluster_tstatsall<= .025);
        else
            display('No significant negative clusters!')
            neg_sig=[];
        end
        if isempty(neg_sig)==0
            neg = ismember(tstats_all.negclusterslabelmat,neg_sig);
        elseif isempty(neg_sig)==1
            
            display('No significant negative clusters!')
        end
    elseif isneg==0
        
        display('No negative clusters!');
        neg_sig=[];
    end
    
    
    %% Plot the Results of the Non-Parametric Cluster-based Test
    if Time_Interval(2)-Time_Interval(1)==timestep
        
        figure;
        
        D = diff(tstats_all.time);     %resolution is generally in the order of 2ms.
        TimeInt_diff=tstats_all.time(end)-tstats_all.time(1);
        SR=DataIn_pp{1,1,1}.fsample;
        timestep1=timestep/(1/SR);
        sample_count=length(tstats_all.time);
        if (tstats_all.time(1)+timestep) > tstats_all.time(end)
            timestep=tstats_all.time(end)-tstats_all.time(1);
            j=[tstats_all.time(1):timestep:tstats_all.time(end)];    %[DataIn{1,1}.time(1):timestep:tstats_all.time(end)]
        else
            
            j=[tstats_all.time(1):timestep:tstats_all.time(end)];
        end
        m=[1:timestep1:sample_count];
        
        
        subplot(1,3,1)
        
        for k1=1:length(m)-1
            
            cfg=[];
            cfg.xlim=[j(k1) j(k1+1)];
            cfg.zlim=[-4 4];
            cfg.marker='on';
            cfg.comment='xlim';
            cfg.commentpos='middletop';
            cfg.layout=layout1;
            ft_topoplotER(cfg,avgSTD);
            colormap(gca,jet);
        end
        
        
        k1=[];
        subplot(1,3,2)
        
        for k1=1:length(m)
            
            cfg=[];
            cfg.xlim=[j(k1) j(k1+1)];
            cfg.zlim=[-4 4];
            cfg.marker='on';
            cfg.markercolor = [1 1 1];
            cfg.markersize = 18;
            cfg.highlight = 'on';
            cfg.comment='xlim';
            cfg.commentpos='middletop';
            %cfg.colorbar='WestOutSide';
            cfg.layout=layout1;
            ft_topoplotER(cfg,avgDEV);
            colormap(gca,jet);
        end
        
        subplot(1,3,3)
        for k=1:length(m)
            
            %subplot(ceil(length(m)/4),4,k)
            cfg=[];
            cfg.xlim=[j(k) j(k+1)];
            
            if isempty(pos_sig)==0 && isempty(neg_sig)==0
                disp('**************Some significant stuff: positive and negative*********************************');
                pos_int = all(pos(:,m(k):m(k+1)),2);
                neg_int = all(neg(:,m(k):m(k+1)),2);
                cfg.highlightchannel = find(pos_int | neg_int);
            elseif isempty(pos_sig)==0 && isempty(neg_sig)==1
                display('No significant negative clusters!');
                pos_int = all(pos(:,m(k):m(k+1)),2);
                cfg.highlightchannel = find(pos_int);
            elseif isempty(pos_sig)==1 && isempty(neg_sig)==0
                display('No significant positive clusters!');
                neg_int = all(neg(:,m(k):m(k+1)),2);
                cfg.highlightchannel=find(neg_int);
            elseif isempty(pos_sig)==1 && isempty(neg_sig)==1
                cfg.highlightchannel=[];
                display('*****************No significant anything!*************************');
            end
            
            
            if ~isempty(find(tstats_all.mask))
                diff_DEVvsSTD_all.avg=diff_DEVvsSTD_all.avg.*tstats_all.mask;  %.*tstats_all.mask;
                cfg.highlight='on';
                cfg.highlightchannel=find(tstats_all.mask);
                cfg.markerfontsize=11;
                cfg.style='both';
                cfg.comment='xlim';
                cfg.commentpos='middletop';
                cfg.layout = layout1;
                cfg.zlim=[-2 2];
                cfg.maskparameter=tstats_all.mask ;
                cfg.maskstyle='saturate';
                ft_topoplotER(cfg,diff_DEVvsSTD_all)
                colormap(gca,jet);
            else
                diff_DEVvsSTD_all.avg = diff_DEVvsSTD_all.avg;  %.*tstats_all.mask;
                cfg.highlight='on';
                cfg.markerfontsize=11;
                cfg.style='both';
                cfg.comment='xlim';
                cfg.commentpos='middletop';
                cfg.layout = layout1;
                cfg.zlim=[-2 2];
                cfg.maskparameter=tstats_all.mask ;
                cfg.maskstyle='saturate';
                ft_topoplotER(cfg,diff_DEVvsSTD_all)
                colormap(gca,jet);
                
                
            end
        end
        
    end
    
    %% Present cluster-based permutation results in time steps.
    if strcmp(timestep_anal,'yes')
        
        scrsz = get(groot,'ScreenSize');
        scrsz(1,4) = scrsz(1,4)/2;
        scrsz(1,3) = scrsz(1,3)/2;
        f2 =figure('Position',scrsz); set(f2,'Color',[1 1 1]);
        %f2 =figure; set(f2,'PaperUnits','normalized','PaperPosition',[0.0235308 0.0272775 0.894169 0.909249],'Color',[1 1 1]);
        
        
        difftime = diff(diff_DEVvsSTD_all.time);    %resolution is generally in the order of 2ms.
        SR = DataIn_pp{1,1,1}.fsample;
        timestep1 = floor(timestep/(1/SR))-1;
        SR = DataIn_pp{1,1,1}.fsample;
        j1= diff_DEVvsSTD_all.time(1):(difftime(2)*timestep1):diff_DEVvsSTD_all.time(end);
        m1 = 1:timestep1:length(diff_DEVvsSTD_all.time);
        
        
        for k1=1:length(m1)-1  %-1
            
            subplot(3,length(m1),k1)  
            cfgds = [];
            cfgds.xlim = [j1(k1) j1(k1+1)];
            cfgds.marker = 'on';
            cfgds.highlight='on';
            cfgds.comment='xlim';
            cfgds.commentpos='middletop';
            cfgds.layout = layout1;
            cfgds.highlightcolor=[0 0 0];
            cfgds.zlim=[-5 5];
            cfgds.shading='interp';
            cfgds.style = 'fill';
            cfgds.colormap=jet;
            
            ft_topoplotER(cfgds,avgDEV);
            
            subplot(3,length(m1),k1+length(m1)) 
            ft_topoplotER(cfgds,avgSTD);
            
            
            subplot(3,length(m1),k1+(length(m1)*2))
            cfg=[];
            cfg.xlim = [j1(k1) j1(k1+1)];
            cfg.zlim = [-2 2];
            
            if isempty(pos_sig)==0 && isempty(neg_sig)==0
                disp('**************Some significant stuff: positive and negative*********************************');
                pos_int = all(pos(:,m1(k1):m1(k1+1)),2);
                neg_int = all(neg(:,m1(k1):m1(k1+1)),2);
                cfg.highlightchannel=find(pos_int | neg_int);
            elseif isempty(pos_sig)==0 && isempty(neg_sig)==1
                display('No significant negative clusters!');
                pos_int=all(pos(:,m1(k1):m1(k1+1)),2);
                cfg.highlightchannel=find(pos_int);
            elseif isempty(pos_sig)==1 && isempty(neg_sig)==0
                display('No significant positive clusters!');
                neg_int = any(neg(:,m1(k1):m1(k1+1)),2);
                cfg.highlightchannel=find(neg_int);
            elseif isempty(pos_sig)==1 && isempty(neg_sig)==1
                cfg.highlightchannel=[];
                display('*****************No significant anything!*************************');
            end
            
            sig_int = zeros(64,1);
            for chancnt = 1:64
                m = find(squeeze(tstats_all.prob(chancnt,:))== min(tstats_all.prob)); %Will return 1 if all time points in the current window are significant (all ones);
                if (length(m)/length(squeeze(tstats_all.prob(chancnt,:)))) <= 0.8
                    sig_int(chancnt) = 0;
                else
                    sig_int(chancnt) = 1;
                end
            end
%       
%             
            assignin('base','sig_int',sig_int);
            
            diff_DEVvsSTD_all.avg = diff_DEVvsSTD_all.avg;   %.*tstats_all.prob;
            diff_DEVvsSTD_all.logprob = log10(tstats_all.prob).*-1;
            cfg.parameter='logprob';
            cfg.marker = 'off';
            cfg.highlight='on';
            cfg.comment='xlim';
            cfg.commentpos='middletop';
            cfg.layout = layout1;
            cfg.highlightsymbol='.';
            cfg.highlightcolor=[0 0 0];
            cfg.zlim=[-5 5];
            cfg.shading='interp';
            cfg.style = 'fill';
            cfg.colormap=jet;
            
            if sum(sig_int)> 0
                cfg.highlightchannel = find(sig_int);
            else
                cfg.highlight = 'off';
            end
            cfg.highlightsize=18;
            ft_topoplotER(cfg,diff_DEVvsSTD_all);
            
            A=unique(log10(tstats_all.prob).*-1);
            colormap(jet);
            %caxis([-5 5])
            caxis([0 1.3])
            
            
        end
        
    end
    
    
    %% PLOT THE RESULTS OF CLUSTER-BASED NON-PARAMETRIC TEST ON INDIVIDUAL ELECTRODES
    
    %reuse the design defined above.
    
    cfg1=[];
    cfg1.method='montecarlo';
    cfg1.statistic= 'ft_statfun_depsamplesT';
    cfg1.correctm='cluster';
    cfg1.clusteralpha=0.05;
    cfg1.clustertail=0;
    cfg1.tail=0;  % two-sided test
    cfg1.correcttail='prob';
    cfg1.alpha=0.025;
    cfg1.numrandomization=1000;
    cfg1.design=design;
    cfg1.ivar=2;    %independent variable
    cfg1.uvar=1;   %unit of observation
    cfg1.channel={DataIn{1,1}.label{1:64}};
    cfg1.neighbours=[];
    
    tstats_single=ft_timelockstatistics(cfg1,DEV{:},STD{:});
    
    f1=figure; set(f1,'Color',[1 1 1]);
    s=zeros(length(cfg1.channel),1);
    st=zeros(length(cfg1.channel),1);
    
    for chan_cnt=1:length(cfg1.channel)
        
        pcurr=tstats_single.prob(chan_cnt,:);
        i=find([pcurr>cfg1.alpha]);
        i2=find([pcurr<=cfg1.alpha]);
        pcurr(i)=0;
        pcurr(i2)=10;
        s(chan_cnt)=subplot(8,8,chan_cnt);
        st(chan_cnt)=stem(s(chan_cnt),tstats_single.time,pcurr);
        set(st(chan_cnt),'MarkerSize',2,'Color','g');
        
        grid on
        xlabel('time (seconds)');
        ylim([0 1])
        xlim([-0.2 1])
        title(s(chan_cnt),DataIn{1,1}.label{chan_cnt});
        
        set(s(chan_cnt),'HitTest','on','SelectionHighlight','on','UserData',{avgSTD.avg(chan_cnt,:),avgDEV.avg(chan_cnt,:),avgSTD.time,tstats_single.time,pcurr,STD_char,DEV_char,DataIn{1,1}.label{chan_cnt},i2},'NextPlot','add');
        set(s(chan_cnt),'ButtonDownFcn',@plotsingle1)
        disp(get(s(chan_cnt),'UserData'));
        
    end
    
    %% CARRY OUT CLUSTER-BASED PERMUTATION TEST ON INDIVIDUAL ELECTRODE DATA
elseif strcmp(anal_type,'single_elec')
    
    statsout=cell(1,channum);
    raweffect=cell(1,channum);
    f1=figure; set(f1,'Color',[1 1 1],'Position',[100 100 1400 1000]);
    s=zeros(channum,1);
    DataIn=squeeze(DataIn);
    
    
    
    for ecnt=1:channum
        
        cfg=[];
        cfg.method='montecarlo';
        cfg.statistic= 'ft_statfun_depsamplesT';  %unit of observation here is subject or time point (if time point means that need to use indepsamplesT)
        cfg.correctm='cluster';
        cfg.clusteralpha=0.05;
        cfg.clustertail=0;
        cfg.clusterstatistic='maxsum';
        cfg.tail=0; % two-sided test
        cfg.correcttail='prob';
        cfg.alpha=0.025;
        cfg.numrandomization=1000;
        cfg.latency=Time_Interval;
        
        design=zeros(2,numsuj*2);  %initialise the design matrix
        design(1,1:numsuj)=1:numsuj;
        design(1,numsuj+1:numsuj*2)=1:numsuj;
        design(2,1:numsuj)=1;
        design(2,numsuj+1:numsuj*2)=2;
        
        cfg.design=design;
        cfg.ivar=2;
        cfg.uvar=1;
        cfg.neighbours=[];
        chan={DataIn{1,1}.label{1:channum}};
        chanind=find(strcmp(chan{ecnt},STD{1,1}.label));
        cfg.channel=chan{chanind};  %define the current electrode of interest
        
        statsout{1,ecnt}=ft_timelockstatistics(cfg,DataIn{:,1},DataIn{:,2});  %DataIn{:,1}= standard     DataIn{:,2}= deviation
        
    end
    assignin('base','statsout', statsout);
    
    indx=find(DataIn{1,1}.time>=statsout{1,1}.time(1) & DataIn{1,1}.time<=statsout{1,1}.time(end));
    Data_elec=cell(length(chan),size(DataIn,2));
    
    for ccnt=1:size(DataIn,2)
        for ecnt=1:length(chan)
            ecnt
            D=zeros(size(DataIn,1),size(DataIn{1,1}.avg,2));
            for scnt=1:size(DataIn,1)
                
                D(scnt,:)=DataIn{scnt,ccnt}.avg(ecnt,:);
                
            end
            Data_elec{ecnt,ccnt}=D;
        end
    end
    
    
    colrs=[ones(1,3).*[0.2 0.7 0.2];ones(1,3).*[0.2 0.2 0.5];ones(1,3).*[0.5 0.2 0.2]];
    CREx_elecconfig_plot(DataIn_eeglab{1,1}(1).chanlocs, chan, Data_elec,statsout{1,1}.time, colrs,statsout,DataIn{1,1}.time,indx,condnames);
    
    
    
end








