function EEGdata = loadEEGdata(RawFilesPath, RawFilesName, StandChanLocsPath, CustomChanLocsPath)
    %% SETTING UP EEG AND PREPARING FOR EEGLAB RELATED FILES
%     clear;
    % Running EEGLAB
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    % Loading a montage (placement of 14 channels of EEG)
    %{  To create a customized montage: typing this in prompt
    %    (Modify the labels according to the experiment)
    %    >> chanlocs = struct('labels', { 'cz' 'c3' 'c4' 'pz' 'p3' 'p4' 'fz' 'f3' 'f4'});
    %    >> pop_chanedit( chanlocs );
    %    Click Ok > Save.ced > Ok
    %   https://sccn.ucsd.edu/wiki/A03:_Importing_Channel_Locations -(Retrieving standardized channel locations) 
    %}

    %STEP 0: Change the option to use double precision
    pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1,...
            'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0,...
            'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1,...
            'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
            'option_checkversion', 1, 'option_chat', 0);
    %% DEFINING VARIABLES RELATED TO FILE

    % Defining a variable that has a path of raw_files
    raw_files = RawFilesPath;

    % Putting all raw files into variable SetFiles (structure array)
    % SetFiles = dir(['Hyper1a_*.hdf5']);
    raw_files_templ_name = RawFilesName; 
    SetFiles = dir([raw_files_templ_name]);

    %% SORTING THE NAME CORRECTLY IN THE STRUCTURE (SETFILES.NAME)

    % extract the numbers
    filenames = {SetFiles.name};   
    filenum = cellfun(@(x)sscanf(x,raw_files_templ_name), filenames);

    % sort them, and get the sorting order
    [~,Sidx] = sort(filenum); 

    % use to this sorting order to sort the filenames
    SetFiles = SetFiles(Sidx);
    close;

    %% Loading files in one batch
    ALLEEG = [];
    for SubjID = 1: length(SetFiles) 
        disp(['Uploading EEG data ' num2str(SubjID) '.' ])
        loadName = SetFiles(SubjID).name;    
        % Loading raw files to EEGLAB
        EEG = pop_loadhdf5('filename',loadName,'filepath',raw_files,'rejectchans',[],'ref_ch',[]);% no pop-up window
        % Update one by one
        [ALLEEG,EEG,CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
        % Loading a standard EEG channel location
        % StandChanLocs = 'C:\\Program Files\\eeglab14_1_2b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp';
        StandChanLocs = StandChanLocsPath;
        % Loading a customized EEG channel locations (16 channels)
        % CustomChanLocs    = 'G:\My Drive\PhD_related_stuff\Codes\Hyperscanning-analysis\channel_location_16.ced';
        CustomChanLocs    = CustomChanLocsPath;
        % Integrating the customized (16 channel) locations into each EEG dataset
        EEG = pop_chanedit(EEG, 'lookup',StandChanLocs,'load',{CustomChanLocs 'filetype' 'autodetect'});
        [ALLEEG,EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end
    EEGdata = ALLEEG;
%     eeglab redraw
end

