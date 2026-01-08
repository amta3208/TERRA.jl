function prob_setup = readProbSetupTERRA(inputFolderPath)
% Read the prob_setup.inp file from TERRA

    % Save current path for later return
    runPath = pwd;

    % Initialize comment characters
    commChars = commCharsTERRA;

    % Variables that are string type
    stringVars = {'DATABASE_PATH','CHEM_FILE_NAME','MOLDATA_PATH'};

    % Open the prob_setup file and load the variables into a structure
    cd(inputFolderPath)
    fullInputFolderPath = pwd;
    fh = fopen('prob_setup.inp');
    while ~feof(fh)
        line = fgetl(fh);
        if (length(strtrim(line)) >= 1)  && ~contains(commChars,line(1))
            splitLine = split(line,'=');
            varName = splitLine{1};
            if any(strcmp(varName,stringVars))
                varVal = strtrim(splitLine{2});         % save as string if var value is string
            else
                varVal = str2double(splitLine{2});      % convert to double if var isn't a string
            end
            prob_setup.(varName) = varVal;              % populate the relevant field in struct
        end
    end

    % Handle the case where databasePath is a relative path
    if ~isfolder(prob_setup.DATABASE_PATH)
      testPath = fullfile(strrep(fullInputFolderPath, '/input', ''), prob_setup.DATABASE_PATH);
      if isfolder(testPath)
        prob_setup.DATABASE_PATH = testPath;
      else
        error('Unable to find database path at: \n%s \n -or- \n%s \n', prob_setup.DATABASE_PATH, testPath)
      end
    end

    % Close the file and return to original path
    fclose(fh);
    cd(runPath)
end
