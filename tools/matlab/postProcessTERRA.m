function caseData = postProcessTERRA(casePath, caseSaveName, debug, overwrite)
% Postprocess a TERRA simulation. Will save the output to a .mat file and,
% if the user requests an output, the postprocessed simulation results will
% be returned as well.

  %% Set debug and overwrite variable if unset
  if (nargin < 3)
    debug = true;
  end
  if (nargin < 4)
    overwrite = false;
  end

  %% Make paths
  inputPath = sprintf('%s/input',casePath);
  outputPath = sprintf('%s/output',casePath);
  statesPath = sprintf('%s/output/states',casePath);
  sourcePath = sprintf('%s/output/sources',casePath);

  %% Read setup data

  fprintf('Reading prob_setup.inp... \n');
  prob_setup = readProbSetupTERRA(inputPath);
  databasePath = prob_setup.DATABASE_PATH;

  fprintf('Reading species.dat... \n');
  species = readSpeciesTERRA(databasePath);

  %% Read output data

  % Initializations
  units = struct();
  labels = struct();

  fprintf('Reading excited state data...\n');
  [levels,species] = readLevelsTERRA(statesPath,species);

  % Read results
  fprintf('Reading result files...\n');
  [result, units, labels] = readResultTERRA(units,labels,outputPath,species,levels,prob_setup);
  result = trimStructTERRA(result,debug);    % Remove extra from the end of files to make them all the same length

  % Read sources
  fprintf('Reading source files...\n');
  [source, labels, units] = readSourceTERRA(units,labels,sourcePath,species,levels,prob_setup);
  source = trimStructTERRA(source,debug);    % Remove extra from the end of files to make them all the same length

  %% Save run file
  % Save the output to a .mat file in the case directory
  saveOutputsTERRA(casePath, caseSaveName, result, source, units, labels, levels, species, prob_setup, overwrite);

  % Assign output if one is requested
  if (nargout == 1)
    caseData = struct('result',result,'source',source,'units',units,'labels',labels,'species',species,'prob_setup',prob_setup);
  end
end
