%%UPDATE 29/04/19
% Gamma di ganti freq 40
%% Isilah titik titik dibawah ini
run E:\eeglab14_1_2b\eeglab; %Masukin path EEGLab
nchannels = 14 ; % Masukan jumlah channel yang ingin dihitung
ChanLocs = 'E:\\eeglab14_1_2b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp';
ChanLocs2 = 'E:\eeglab14_1_2b\emotivupdated.ced';
ChanNum = [3:16]; %if Emotiv
% ChanNum = [1:19]; %if Deymed / 19 Channels Stuffs
JumlahDomain = 8; %Obvious lah
CutPath = 'C:\Users\DAI 02 - Neurolab\Desktop\Cutting Folder\CutFiles\'; %Masukin Path file yang udah di cut
close;
%% Yang ini tidak usah di isi
%Cut Frame on 1366x768
deltaFrame = [181 455 154 154];
thetaFrame = [401 455 154 154];
alphaFrame = [622 455 154 154];
betaFrame  = [842 455 154 154];
gammaFrame = [1062 455 154 154];

cd(CutPath);
mkdir PP;
ProcPath = [CutPath 'PP\'];
cd(ProcPath);
mkdir BaseFile;
mkdir Plot;
ProcBasePath = [ProcPath 'BaseFile\'];
ImgPath = [ProcPath 'Plot\'];

SetFiles = dir([CutPath '*.set']);
ProcFiles = dir([ProcPath '*.set']);
ProcFDT = dir([ProcPath '*.fdt']);
EdfFiles = dir([CutPath '*.edf']);


%STEP 0: Change the option to use double precision
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1,...
        'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0,...
        'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1,...
        'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
        'option_checkversion', 1, 'option_chat', 0);

%% PREPROCESSING
for SubjID = 1: length(SetFiles);
    loadName = SetFiles(SubjID).name;
    dataName = loadName(1:end-4);
    saveName = [dataName '_PP'] %adding 'PP' (PreProcessed) to the end 
    
    %Step 1: Import data. REVISED 
  % EEG = pop_biosig(loadName, CutPath); %IF .EDF
    EEG = pop_loadset(loadName, CutPath); %IF .SET
    EEG.setname = dataName;

    %STEP 2: Changing sample rate
    EEG = pop_resample(EEG, 250);

    %STEP 3: High-pass filter and Low pass filter
    EEG = pop_eegfiltnew(EEG,1,45,826);

    %STEP 4: Select channels
    EEG = pop_select( EEG,'channel',ChanNum);
    
    %STEP 5: Importing channel location
    EEG = pop_chanedit(EEG, 'lookup',ChanLocs,'load',{ChanLocs2 'filetype' 'autodetect'});
    
        %Keeping original EEG file
        originalEEG = EEG; % Just back up original EEG data
        
    %STEP 6: Removing bad channels by using raw_clean plug in
    EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5);

    %STEP 7: Interpolate channels.
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');

    %STEP 8: Apply average reference after adding initial reference
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    EEG = pop_select( EEG,'nochannel',{'initialReference'});

    %STEP 9: Epoching data 1 to 3 sec
    EEG = eeg_regepochs(EEG, 'limits', [1 3] , 'extractepochs', 'on'); % RUBAH timing epoch

    %STEP 10: Automatic epoch rejection
    EEG = pop_autorej(EEG, 'threshold', 1000,'startprob',5,'maxrej', 5, 'nogui','on'); %nogui diganti on, eegplot hapus aja
    
    %STEP 11: Rejection epoch by probability (6SD single channel, 2SD for all channels)
    EEG = eeg_checkset( EEG );
    EEG = pop_jointprob(EEG,1,[1:14] ,6,2,1,0,0,[],0);
    
    %STEP 12: Running ICA
    EEG = eeg_checkset( EEG );
    EEG = pop_runica(EEG, 'extended',1,'interupt','on');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    %STEP 13 : Checking whether EEG data contains ICA decomposition
    EEG = eeg_checkset(EEG, 'ica');
    
    %STEP 14: Removing line noise with cleanLine;
    EEG = pop_cleanline(EEG,'Bandwidth',2,'ChanCompIndices',[1:EEG.nbchan],...
        'ComputeSpectralPower',0,'LineFrequencies',[60 120],...
        'NormalizeSpectrum',0,'LineAlpha',0.01,'PaddingFactor',2,...
        'PlotFigures',0,'ScanForLines',1,'SignalType','Channels',...
        'SmoothingFactor',100,'VerboseOutput',1); % HAPUS SlidingWinLength dan SlidingWinstep
    
    %STEP 15: Rejecting ICA by extreme value
    EEG = pop_eegthresh(EEG,0,[1:5] ,-20,20,1,2.996,0,1);

    %STEP 16: Finding Power for each frequency
    %%%%%%%%%%%%%%%%%%% Finding power for each frequency band %%%%%%%%%%%%%%%%%%%
    
    EEG = pop_saveset(EEG, 'filename',saveName,'filepath', ProcPath);  
end
%% Statistic
for ProcID = 1:length(ProcFiles)
    loadProc = ProcFiles(ProcID).name;
    procData = loadProc(1:end-4);
    
    EEG = pop_loadset(loadProc, ProcPath); %IF .SET
    EEG.setname = procData;
    
    for n = 1: nchannels

        [spectra,freqs] = spectopo(EEG.data(n,:,:), 0, EEG.srate, 'plot', 'off'); % Sesuaikan channel mana yang mau diambil

        % delta=1-4, theta=4-8, alpha=8-13, beta=13-30, gamma=30-80
        deltaIdx{n} = find(freqs>1 & freqs<=4);
        thetaIdx{n} = find(freqs>4 & freqs<=8);
        alphaIdx{n} = find(freqs>8 & freqs<=13);
        betaIdx{n}  = find(freqs>13 & freqs<=30);
        gammaIdx{n} = find(freqs>30 & freqs<=80);

        % compute absolute power
        deltaPower{n} = mean(10.^(spectra(deltaIdx{n})/10));
        thetaPower{n} = mean(10.^(spectra(thetaIdx{n})/10));
        alphaPower{n} = mean(10.^(spectra(alphaIdx{n})/10));
        betaPower{n}  = mean(10.^(spectra(betaIdx{n})/10));
        gammaPower{n} = mean(10.^(spectra(gammaIdx{n})/10));
    end
    
        deltaCell= cell2mat(deltaPower);
        thetaCell= cell2mat(thetaPower);
        alphaCell= cell2mat(alphaPower);
        betaCell = cell2mat(betaPower);
        gammaCell= cell2mat(gammaPower);
        
        % Dikumpul per wave
        deltaAll(ProcID,:) = deltaCell
        thetaAll(ProcID,:) = thetaCell
        alphaAll(ProcID,:) = alphaCell
        betaAll (ProcID,:) = betaCell
        gammaAll(ProcID,:) = gammaCell
        
end
%Indexing p-value per Domain
for StatID = 1:JumlahDomain+1:length(ProcFiles);
%     loadStat = ProcFiles(StatID).name;
%     statData = loadStat(1:end-4);
     for Stat = 1:JumlahDomain;
     [h,p] = ttest(deltaAll(StatID,:),deltaAll(StatID+Stat,:));
     PValDelta(StatID+Stat) = p ;
     end
     for Stat = 1:JumlahDomain;
     [h,p] = ttest(thetaAll(StatID,:),thetaAll(StatID+Stat,:));
     PValTheta(StatID+Stat) = p;
     end
     for Stat = 1:JumlahDomain;
     [h,p] = ttest(alphaAll(StatID,:),alphaAll(StatID+Stat,:));
     PValAlpha(StatID+Stat) = p;
     end
     for Stat = 1:JumlahDomain;
     [h,p] = ttest(betaAll(StatID,:),betaAll(StatID+Stat,:));
     PValBeta(StatID+Stat) = p;
     end
     for Stat = 1:JumlahDomain;
     [h,p] = ttest(gammaAll(StatID,:),gammaAll(StatID+Stat,:));
     PValGamma(StatID+Stat) = p;
     end
end
     %removing zeroes for Indexing
PValDelta(PValDelta==0) =[];
PValTheta(PValTheta==0) =[];
PValAlpha(PValAlpha==0) =[];
PValBeta(PValBeta==0) =[];
PValGamma(PValGamma==0) =[];

%Dikumpul jadi satu workspace also for 1 row Indexing
PValSubj = [PValDelta; PValTheta; PValAlpha; PValBeta; PValGamma]
[PValMin,PValIdx] = min(PValSubj)

%% Plotting
for z = 1:length(PValIdx);
    loadStat = ProcFiles(z).name;
    statData = loadStat(1:end-4);
 
%indexing smallest p-value sebagai domain (1=delta - 5=gamma)
 if PValIdx(z) == 1
           Coord(z) = {deltaFrame}
           ImgName{z} = ['_delta.png']
        elseif PValIdx(z) == 2
           Coord(z) =  {thetaFrame}
           ImgName{z} = ['_theta.png']
        elseif PValIdx(z) == 3
           Coord(z) =  {alphaFrame}
           ImgName{z} = ['_alpha.png']
        elseif PValIdx(z) == 4
          Coord(z) =  { betaFrame}
          ImgName{z} = ['_beta.png']
  elseif PValIdx(z) == 5
          Coord(z) = {gammaFrame}
          ImgName{z} = ['_gamma.png']
          
 end
end
%move BaseFile to new Folder
for BaseSet = 1:JumlahDomain+1:length(ProcFiles)
movefile(ProcFiles(BaseSet).name, ProcBasePath)
movefile (ProcFDT(BaseSet).name, ProcBasePath)
end
ProcFiles = dir([ProcPath '*.set']); %refresh filelist in the folder after moved
  for PlotID = 1:length(ProcFiles)
      PlotFile = ProcFiles(PlotID).name
      PlotName = [PlotFile(1:end-4),ImgName{PlotID}]
            EEG = pop_loadset(PlotFile, ProcPath);
      EEG.setname = PlotName;
        %open figure&
        figure('units','normalized','outerposition',[0 0 1 1]);  %%Full Screen
        title(PlotName); %Title
        % Plotting Head
        pop_spectopo(EEG, 1, [], 'EEG', ...
        'percent', 100, 'freq', [2 6 12 25 40], 'freqrange',[0 45], 'electrodes','off', 'overlap', 0);
        set(gcf,'color','w');
        %Theta
        F = getframe(gcf, Coord{PlotID});
        figure;
        imshow(F.cdata);
        Image = frame2im(F);
               imwrite(Image,[ImgPath PlotName],'png')
        close
        close
   end