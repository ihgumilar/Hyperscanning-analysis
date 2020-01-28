% Objective: Chunking the data into pre-training(trigger 1) and post-training(trigger 3)
% Progress : 1. Loading original files (*.hdf5)
%            2. Cut file into 120 seconds duration (pre-training) 
%            3. Cut file into 120 seconds duration (post-training)
%            4. Pre-processing files, 1st pair not working (finding why ?)
%
% Future development : 1. Filter the data with different frequency resolution for each band
%% Step 1 : Load original raw files of EEG (*.hdf5)
clear;
% Running EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Change the option to use double precision
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1,...
        'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0,...
        'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1,...
        'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
        'option_checkversion', 1, 'option_chat', 0);

% Defining a variable that has a path of raw_files
cd 'G:\My Drive\PhD_related_stuff\Codes\Hyperscanning-analysis'           % CHANGE HERE
raw_files = 'G:\My Drive\PhD_related_stuff\Codes\Hyperscanning-analysis'; % CHANGE HERE

% Putting all raw files into variable SetFiles (structure array)
SetFiles = dir(['Hyper1a_*.hdf5']);

% Sorting the name correctly in the structure (Setfiles.name)
% extract the numbers
  filenames = {SetFiles.name};   
  filenum = cellfun(@(x)sscanf(x,'Hyper1a_%d.hdf5'), filenames);

% sort them, and get the sorting order
  [~,Sidx] = sort(filenum); 

% use to this sorting order to sort the filenames
  SetFiles = SetFiles(Sidx);
  close;

% Loading the raw files (*.hdf5) in one batch
ALLEEG = [];
for SubjID = 1: length(SetFiles) 
    loadName = SetFiles(SubjID).name;    
    % Loading raw files to EEGLAB
    EEG = pop_loadhdf5('filename',loadName,'filepath',raw_files,'rejectchans',[],'ref_ch',[]);% no pop-up window
    % Update one by one
    [ALLEEG,EEG,CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
    % Loading a standard EEG channel location
    StandChanLocs = 'C:\\Program Files\\eeglab14_1_2b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp';
    % Loading a customized EEG channel locations (16 channels)
    CustomChanLocs    = 'G:\My Drive\PhD_related_stuff\Codes\Hyperscanning-analysis\channel_location_16.ced';
    % Integrating the customized (16 channel) locations into each EEG dataset
    EEG = pop_chanedit(EEG, 'lookup',StandChanLocs,'load',{CustomChanLocs 'filetype' 'autodetect'});
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end
% eeglab redraw
%% Step 2 : Cut only the 120 seconds of pre-training files (trigger 1) & save them as *.set

for SubjID = 1: length(SetFiles) 
    % Load the raw file one by one
    EEG = ALLEEG(SubjID);
    % Cut the file from trigger 1 (120 seconds)
    EEG = pop_rmdat( EEG, {'Trigger 1'},[1 120] ,0);
    % Set the name for files that have been cut off
    saveName_pre = [EEG.setname '_Pre']; % CHANGE HERE
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, SubjID,'setname',saveName_pre,'gui','off'); 
    % Save the files
    EEG = pop_saveset( EEG, 'filename',saveName_pre,'filepath','G:\\My Drive\\PhD_related_stuff\\Codes\\Exp1_preTraining'); % CHANGE HERE
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);   
end
% eeglab redraw
%% Step 3 : Cut only the 120 seconds of post-training files (trigger 3) & save them as *.set

for SubjID = 1: length(SetFiles) 
    % Load the raw file one by one
    EEG = ALLEEG(SubjID);
    % Cut the file from trigger 3 (120 seconds)
    EEG = pop_rmdat( EEG, {'Trigger 3'},[1 120] ,0);
    % Set the name for files that have been cut off
    saveName_post = [EEG.setname '_Post']; % CHANGE HERE
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, SubjID,'setname',saveName_post,'gui','off');
    % Save the files
    EEG = pop_saveset( EEG, 'filename',saveName_post,'filepath','G:\My Drive\PhD_related_stuff\Codes\Exp1_posTraining'); % CHANGE HERE
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);   
end
ALLEEG(1:24) = []; % Delete the original raw files (*.hdf5)
eeglab redraw

%% Pre-processing the data.
% (NOTE !) The 1st pair was excluded. Try to run with the whole data and find the issue with the first pair 

for i = 1:length (ALLEEG) % We can adjust this according to the number of our data
    disp(['Processing file...iteration ', num2str(i)])
   % EEG1 is the first dataset of the pair, whilst EEG2 is the second
    [EEG, ALLEEG, CURRENTSET] = eeg_retrieve(ALLEEG,i*2-1);
    saveName1 = [EEG.setname '_PP']; %adding 'PP' (PreProcessed) to the end 
    EEG1 = EEG;

    [EEG, ALLEEG, CURRENTSET] = eeg_retrieve(ALLEEG,i*2);
    saveName2 = [EEG.setname '_PP']; %adding 'PP' (PreProcessed) to the end 
    EEG2 = EEG;

    % STEP 1: High-pass filter and Low pass filter both
    EEG1 = pop_eegfiltnew(EEG1,1,[],1690);
    EEG2 = pop_eegfiltnew(EEG2,1,[],1690);


    % Keeping original EEG files
    originalEEG1 = EEG1;
    originalEEG2 = EEG2;


    % STEP 2: Removing line noise with cleanLine;
    % electricity interference etc (requires plugin)
    EEG1 = pop_cleanline(EEG1,'Bandwidth',2,'ChanCompIndices',[1:EEG.nbchan],...
        'ComputeSpectralPower',0,'LineFrequencies',[60 120],...
        'NormalizeSpectrum',0,'LineAlpha',0.01,'PaddingFactor',2,...
        'PlotFigures',0,'ScanForLines',1,'SignalType','Channels',...
        'SmoothingFactor',100,'VerboseOutput',1); 

    EEG2 = pop_cleanline(EEG2,'Bandwidth',2,'ChanCompIndices',[1:EEG.nbchan],...
        'ComputeSpectralPower',0,'LineFrequencies',[60 120],...
        'NormalizeSpectrum',0,'LineAlpha',0.01,'PaddingFactor',2,...
        'PlotFigures',0,'ScanForLines',1,'SignalType','Channels',...
        'SmoothingFactor',100,'VerboseOutput',1); 


    % STEP 3: clean_rawdata (requires plugin)

    % Keep copies for later
    EEG1_a = EEG1;
    EEG2_a = EEG2;

    % First run of clean_rawdata (work out which points are removed for each dataset)
    % The idea here is that we want to remove the same data in each dataset
    % so that we have clean data on both participants to conduct analysis.
    % (Probably a better way to do this but here is my workaround anyway).
    EEG1 = clean_rawdata(EEG1, 5, [0.25 0.75], 0.8, 4, 20, 0.1);
    EEG2 = clean_rawdata(EEG2, 5, [0.25 0.75], 0.8, 4, 20, 0.1); 
    % Here a fairly aggressive time window rejection parameter of 0.1 is
    % used, and as a result our 120 second datasets are reduced to as low
    % as 30 seconds of really clean data for analysis (after further
    % data rejection steps later on in steps 6 and 7). Might be good to
    % trial different parameters here (depends on the quality of the
    % experiment apparatus, etc), but our results look good so we haven't
    % made changes for now. Need to read clean_rawdata documentation to 
    % fully understand.

    % Combine clean_sample_masks of both datasets (mask of zeros and ones
    % indicating time points that have been rejected/kept)
    mask = EEG1.etc.clean_sample_mask;
    mask(~EEG2.etc.clean_sample_mask)=0;

    % Create set of ranges to discard from both datasets after second run
    % of clean_rawdata (ranges formatted [min1,max1;min2,max2;...])
    compare = 1;
    ranges = [];
    k = 1;
    for j = 1:length(mask)
        if mask(j)==0 & compare==1
            ranges(k,1) = j;
            compare = 0;
        elseif mask(j)==1 & compare==0
            ranges(k,2) = j-1;
            compare = 1;
            k = k+1;
        end
    end

    % Run clean_rawdata again (set time window rejection parameter to
    % 'off' so that no data is rejected)
    EEG1 = clean_rawdata(EEG1_a, 5, [0.25 0.75], 0.8, 4, 20, 'off'); % eeg is dirty we clean it here using raw_clean plugin
    EEG2 = clean_rawdata(EEG2_a, 5, [0.25 0.75], 0.8, 4, 20, 'off'); % eeg is dirty we clean it here using raw_clean plugin

    % Interpolate bad channels using spherical method
    EEG1 = pop_interp(EEG1,originalEEG1.chanlocs,'spherical');
    EEG2 = pop_interp(EEG2,originalEEG2.chanlocs,'spherical');

    % Reject data using ranges matrix (this rejects the same data on both 
    EEG1 = pop_select(EEG1,'nopoint',ranges);
    EEG2 = pop_select(EEG2,'nopoint',ranges);

    % Update clean_sample_mask for both datasets
    EEG1.etc.clean_sample_mask = mask;
    EEG2.etc.clean_sample_mask = mask;

    % STEP 4: Apply average reference after adding initial reference
    EEG1.nbchan = EEG1.nbchan+1;
    EEG1.data(end+1,:) = zeros(1, EEG1.pnts);
    EEG1.chanlocs(1,EEG1.nbchan).labels = 'initialReference';
    EEG1 = pop_reref(EEG1, []); 
    EEG1 = pop_select( EEG1,'nochannel',{'initialReference'});

    EEG2.nbchan = EEG2.nbchan+1;
    EEG2.data(end+1,:) = zeros(1, EEG2.pnts);
    EEG2.chanlocs(1,EEG2.nbchan).labels = 'initialReference';
    EEG2 = pop_reref(EEG2, []); 
    EEG2 = pop_select( EEG2,'nochannel',{'initialReference'});


    % STEP 5: Epoching data into continuous 1 second sections
    EEG1 = eeg_regepochs(EEG1,'limits', [0 1] , 'extractepochs', 'on'); 
    EEG2 = eeg_regepochs(EEG2,'limits', [0 1] , 'extractepochs', 'on');


    % STEP 6: Perform automatic epoch rejection 
    [EEG1 rmepochs] = pop_autorej(EEG1, 'threshold', 500,'startprob',5,'maxrej', 5,'nogui','on');
    EEG2 = pop_select(EEG2,'notrial',rmepochs); % Again, ensure same epochs are rejected in dataset of partner
    [EEG2 rmepochs] = pop_autorej(EEG2, 'threshold', 500,'startprob',5,'maxrej', 5,'nogui','on');
    EEG1 = pop_select(EEG1,'notrial',rmepochs);


    % STEP 7: Further automatic epoch rejection by probability 
    % (6SD single channel, 2SD for all channels)
    EEG1 = eeg_checkset( EEG1 ); 
    [EEG1, locthresh, globthresh, nrej, rej] = pop_jointprob(EEG1,1,[1:16],6,2,1,1,1,[],0);
    rmepochs = find(rej);
    EEG2 = pop_select(EEG2,'notrial',rmepochs);
    EEG2 = eeg_checkset( EEG2 ); 
    [EEG2, locthresh, globthresh, nrej, rej] = pop_jointprob(EEG2,1,[1:16],6,2,1,1,1,[],0);
    rmepochs = find(rej);
    EEG1 = pop_select(EEG1,'notrial',rmepochs);
    
    % Saving files that have been pre-processed
    pathname1 = ['G:\My Drive\PhD_related_stuff\Codes\PreProcessed\' saveName1 '.set']
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG1, 1,'setname',saveName1,'savenew',pathname1,'gui','off');
    pathname2 = ['G:\My Drive\PhD_related_stuff\Codes\PreProcessed\' saveName2 '.set']
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG2, 1,'setname',saveName2,'savenew',pathname2,'gui','off');
    
end
eeglab redraw


%% Save the pre-processed files into text file
% (NOTE !) Make sure leave only the files that have been pre-processed with the
% following code:  ALLEEG(1: n) = [];
% Begin to export into text file
