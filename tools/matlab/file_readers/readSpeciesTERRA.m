function species = readSpeciesTERRA(databasePath)
% Read species.dat file for TERRA from a database folder

    % Save current path for later return
    runPath = pwd;

    % Initialize comment characters
    commChars = commCharsTERRA;

    % Enter the database folder
    cd(databasePath)

    % Read file and extract relevant variables
    fh = fopen('species.inp');
    s = 0;
    while ~feof(fh)
        line = fgetl(fh);
        if (length(line) >= 1)  && ~contains(commChars,line(1))
            s = s + 1;
            splitLine = split(line);
            species.spnm{s}  = strrep(splitLine{1},'+','p');
            if strcmp(species.spnm{s},'E-')
              species.spnm{s} = 'E';
            end
            species.spwt(s)  = str2double(splitLine{2});
            species.ih(s)    = str2double(splitLine{3});
            species.ie(s)    = str2double(splitLine{4});
            species.ies(s)   = str2double(splitLine{5});
            species.h0(s)    = str2double(splitLine{6});
            species.eion(s)  = str2double(splitLine{7});
            species.factr(s) = str2double(splitLine{8});
            species.theta(s) = str2double(splitLine{9});
        end
    end
    fclose(fh);
    cd(runPath)
end
