%% PRE-PROCESSING script
% This script applies pre-processing methods to the EEG data of people 
% participating in the hyperscanned finger-pointing experiment to remove
% sources of artifacts and noise, and reject periods of time too noisey for
% analysis. Before running this script, 48 datasets are loaded into eeglab;
% The first 24 datasets are the raw 120-second pre-training sections of the
% participants, and the second 24 datasets are the raw 120-second
% post-training sections of the participants. Participants 1 and 2 are
% partners in the experiment, as are participants 3 and 4, as are 
% participants ..., as are particpants 23 and 24.

% The pre-processing steps I have used come mainly from Makoto's
% pre-processing pipeline, with modifications to suit our purposes

% Author: Gus Stone

% Loop through pairs of simultaneous pre-training datasets (i=1:12), then
% pairs of simultaneous post-training datasets (i=13:24)
for i = 1:24
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
    
    
    % STEP 8: Running ICA
    % Separates the independent components for each dataset.
    % Afterwards, I decided that it would be best to manually inspect the
    % independent components in order to identify and reject components
    % that were obvious sources of artifacts (eg. ocular movement/eye-blink
    % components), rather than using automatic functions to rejects IC's. 
    % Therefore, the last stage of pre-processing is done manually in 
    % eeglab, using tools -> reject data using ICA -> reject components by 
    % map. See http://mikexcohen.com/lectures.html for help here.
    EEG1 = eeg_checkset( EEG1 );
    EEG1 = pop_runica(EEG1, 'extended',1,'interupt','on'); 

    EEG2 = eeg_checkset( EEG2 );
    EEG2 = pop_runica(EEG2, 'extended',1,'interupt','on'); 

    EEG1 = eeg_checkset(EEG1, 'ica');
    EEG2 = eeg_checkset(EEG2, 'ica');
    
    % Save files (remember to inspect independent components before
    % analysing data)
    pathname1 = ['D:\Github\Hyperscanning\Testing1' saveName1 '.set']
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG1, 1,'setname',saveName1,'savenew',pathname1,'gui','off');
    pathname2 = ['D:\Github\Hyperscanning\Testing1' saveName2 '.set']
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG2, 1,'setname',saveName2,'savenew',pathname2,'gui','off');
end

eeglab redraw