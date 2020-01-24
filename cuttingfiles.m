%% Chunking the data into pre-training(trigger 1) and post-training(trigger 3)
% Step 1 : Load original raw files of EEG

% Step 2 : Cut only the 120 seconds of pre-training files
for SubjID = 1: length(SetFiles) 
    % Calling EEG data one-by-one that is contained in ALLEEG 
    
%     % Cutting files 120-seconds after Trigger 1 (pre-training)
%     EEG = pop_rmdat( EEG, {'Trigger 1'},[1 120] ,0);
%     % Creating a new name for a cuted file
%     saveName_pre = [EEG.setname '_Pre']; 
%     % Seeting path to save
%     savepath = ['G:\My Drive\PhD_related_stuff\Codes' saveName_pre '.set'];
%     % Saving files with new file names
%     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, SubjID,'gui','off');  
    EEG = ALLEEG(SubjID);
    EEG = pop_rmdat( EEG, {'Trigger 1'},[1 120] ,0);
    saveName_pre = [EEG.setname '_Pre'];
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, SubjID,'setname',saveName_pre,'gui','off'); 
    EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath','G:\\My Drive\\PhD_related_stuff\\Codes');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);   
end
eeglab redraw

% Step 3 : Cut only the 120 seconds of post-training files

%% Loading the files into EEGLAB GUI

%% Pre-processing the data (Copy and paste the code here)

%% Save the pre-processed files into text file
