function [rawData, latexVarNames] = readOutputFileTERRA(outputPath,fileName)
% Read an output file from a TERRA simulation and assign it to a single
% double array containing the raw data

  % Save current path for later return
  runPath = pwd;

  % Change directory to the output path
  cd(outputPath);

  % Open the file and skip the top level headers
  fh = fopen(fileName);
  ensureNextLineBeginsWith('TITLE=');
  ensureNextLineBeginsWith('FILETYPE=');

  % Read the variable names
  ensureNextLineBeginsWith('VARIABLES=');
  line = fgetl(fh);
  nvars = 0;
  while ~strcmp2(line,'ZONE T=','start')
    nvars = nvars + 1;
    line = strrep(line,'",','');
    line = strrep(line,'"','');
    tecplotVarName = strtrim(line);
    latexVarNames{nvars} = tecplot2latex(tecplotVarName);
    line = fgetl(fh);
  end
  ensureNextLineBeginsWith('ZONETYPE = ORDERED, DATAPACKING = POINT')

  % Read the data
  rawData = nan(1,nvars);
  id = 0;
  while ~feof(fh)
    if ~feof(fh)
      try
        lineData = readLineData(fh);
      catch
        continue
      end
    end
    id = id + 1;
    lwt = numel(lineData);
    rawData(id,1:lwt) = lineData;
    while lwt < nvars
      try
        lineData = readLineData(fh);
      catch
        lwt = nvars;
      end
      lw = numel(lineData);
      lwt = lwt + lw;
      rawData(id,lwt-lw+1:lwt) = lineData;
    end
  end

  % Return to calling folder
  cd(runPath);

  % Utility functions

  % Ensure header is as expected
  function ensureNextLineBeginsWith(str)
      if ~strcmp2(fgetl(fh),str,'start')
        error('Expected ''%s'' at the beginning of %s.', str, fileName);
      end
  end

  function latexStr = tecplot2latex(tpStr)
    % Subscripts and superscripts
    tpStr = strrep(tpStr,'<sub>','_{');
    tpStr = strrep(tpStr,'</sub>','}');
    tpStr = strrep(tpStr,'<sup>','^{');
    tpStr = strrep(tpStr,'</sup>','}');
    % Greek letters
    tpStr = strrep(tpStr,'<greek>g</greek>','\gamma');
    tpStr = strrep(tpStr,'<greek>r</greek>','\rho');
    tpStr = strrep(tpStr,'<greek>m</greek>','\mu');
    tpStr = strrep(tpStr,'<greek>e</greek>','\epsilon');
    % Math symbols
    tpStr = strrep(tpStr,'<math>a</math>','\langle');
    tpStr = strrep(tpStr,'<math>q</math>','\rangle');
    % Assign output
    latexStr = tpStr;
  end

% Read data line
  function data = readLineData(fh)
    data = str2double(split(strtrim(fgetl(fh))))';
  end
end
