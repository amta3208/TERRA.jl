function [result, units, labels] = readResultTERRA(units,labels,outputPath,species,levels,prob_setup)
% Read result files from a TERRA run

  % Initialize the result structure
  result = struct();

  % Extract useful variables
  spnm = species.spnm;
  mex = species.mex;
  NS = length(mex);
  ies = species.ies;
  ih = species.ih;

  % Create result type fields
  resultNames = {'ener',  'spgm',  'temp', 'flow', 'time'};   % basic types, common to all simulations
  if (prob_setup.ND == 1)
    resultNames{end+1} = 'dist';                              % distance if 1D
  end
  for sp = 1:NS
    if (ies(sp) == 1)
      resultNames{end+1} = sprintf('exgm%d',sp);              % excited state concentrations
      for iex = 1:mex(sp)
        resultType = sprintf('vxgm%d_%d',sp,iex);
        resultFile = sprintf('result-%s.dat',resultType);
        resultPath = makePath(outputPath,resultFile);
        if isfile(resultPath)
          resultNames{end+1} = resultType;
          species.ivs(iex,sp) = 1;
        else
          species.ivs(iex,sp) = 0;
        end
      end
    end
  end
  ivs = species.ivs;
  if (any(ies))
    resultNames{end+1} = 'eion';
  end

  for rf = 1:numel(resultNames)
    resultType = resultNames{rf};

    % Create the result file name
    resultFile = sprintf('result-%s.dat',resultType);
    resultPath = sprintf('%s/%s',outputPath,resultFile);

    % Check if the file exists!
    if ~isfile(resultPath)
      return;
    end

    % Read raw data from output file and reformat into nicely named fields

    switch resultType(1:4)
      case 'time'
        [rawData, latexVarNames] = readOutputFileTERRA(outputPath,resultFile);
        assignStructFields('t',fileData('t [\mus]'),'mu-s','t');
        assignStructFields('dt',fileData('dt [\mus]'),'mu-s','dt');
      case 'dist'
        [rawData, latexVarNames] = readOutputFileTERRA(outputPath,resultFile);
        assignStructFields('x',fileData('x [cm]'),'cm','x');
        assignStructFields('dx',fileData('dx [cm]'),'cm','dx');
      case 'ener'
        [rawData, latexVarNames] = readOutputFileTERRA(outputPath,resultFile);
        assignStructFields('eeex',fileData('e_{eex} [erg/g]'),'erg/g','e_{eex}');
        assignStructFields('erot',fileData('e_{rot} [erg/g]'),'erg/g','e_{rot}');
        assignStructFields('evib',fileData('e_{vib} [erg/g]'),'erg/g','e_{vib}');
        for isp = 1:NS
          if ies(isp)
            fnm = sprintf('e_{ex,%s} [erg/g]',strrep(spnm{isp},'p','+'));
            assignStructFieldsPerSpecies('eex',spnm{isp},fileData(fnm),'erg/g','e_{ex}');
            for iex = 1:mex(isp)
              if ivs(iex,isp)
                fnm = sprintf('e_{v,%s(%d)} [erg/g]',strrep(spnm{isp},'p','+'),iex);
                assignStructFieldsPerSpecies('evx',sprintf('%s_%d',spnm{isp},iex),fileData(fnm),'erg/g','e_{vx}');
              end
            end
          end
        end
      case 'spgm'
        [rawData, latexVarNames] = readOutputFileTERRA(outputPath,resultFile);
        assignStructFieldsAllSpecies('spgam',rawData,'mol/g','\\gamma');
      case 'exgm'
        i = str2double(resultType(5:end));
        [rawData, latexVarNames] = readOutputFileTERRA(outputPath,resultFile);
        assignStructFieldsPerSpecies('exgam',spnm{i},rawData,'mol/g','\\gamma');
      case 'vxgm'
        usid = find(resultType=='_');
        isp = str2double(resultType(5:usid-1));
        iex = str2double(resultType(usid+1:end));
        [rawData, latexVarNames] = readOutputFileTERRA(outputPath,resultFile);
        assignStructFieldsPerSpecies('vibgam',sprintf('%s_%d',spnm{isp},iex),rawData,'mol/g','\\gamma');
      case 'temp'
        [rawData, latexVarNames] = readOutputFileTERRA(outputPath,resultFile);
        assignStructFields('Tt',fileData('T_{t} [K]'),'K','T_t');
        assignStructFields('Teex',fileData('T_{eex} [K]'),'K','T_{eex}');
        assignStructFields('Trot',fileData('T_{rot} [K]'),'K','T_{rot}');
        assignStructFields('Tvib',fileData('T_{vib} [K]'),'K','T_{vib}');
        for isp = 1:NS
          if ies(isp)
            fnm = sprintf('T_{ex,%s} [K]',strrep(spnm{isp},'p','+'));
            assignStructFieldsPerSpecies('Tex',spnm{isp},fileData(fnm),'K','T_{ex}');
            for iex = 1:mex(isp)
              if ivs(iex,isp)
                fnm = sprintf('T_{v,%s(%d)} [K]',strrep(spnm{isp},'p','+'),iex);
                assignStructFieldsPerSpecies('Tvx',sprintf('%s_%d',spnm{isp},iex),fileData(fnm),'K','T_{vx}');
              end
            end
          end
        end
      case 'flow'
        [rawData, latexVarNames] = readOutputFileTERRA(outputPath,resultFile);
        assignStructFields('p',torr2Pa(fileData('p [torr]')),'Pa','p');
        assignStructFields('u',1000*fileData('u [km/s]'),'m/s','u');
        assignStructFields('rho',fileData('\rho [g/cm^{3}]'),'g/cm^{3}','\\rho');
        assignStructFields('avmw',fileData('M_{avg} [g/mol]'),'g/mol','\bar{M}');
      case 'eion'
        [rawData, latexVarNames] = readOutputFileTERRA(outputPath,resultFile);
        for isp = 1:NS
          ispnm_true = strrep(strrep(spnm{isp},'p','+'),'n','-');
          if ~ies(isp) || (ih(isp) == 0)
            continue
          end
          fnm_eps_eex = sprintf('(\\langle\\epsilon_{eex}\\rangle_{d}/I_{0})_{%s} [~]',ispnm_true);
          assignStructFields(sprintf('eps_eex_eion_%s',spnm{isp}),fileData(fnm_eps_eex),'~',strrep(fnm_eps_eex,' [~]',''))
        end
    end

    % Create number density fields and mole fraction fields
    % Delete spgam and exgam afterward - number density is better
    if (isfield(result,'avmw') && isfield(result,'spgam') && ...
       ~isfield(result,'smgam'))
      smgam = 1./(result.avmw);
      ids = 1:min(length(result.spgam.(spnm{1})),length(smgam));
      for i = 1:NS
        spgam = result.spgam.(spnm{i})(ids);
        X = spgam ./ smgam(ids);
        n = spgam .* avnCGS .* result.rho(ids);
        assignStructFieldsPerSpecies('X',spnm{i},X,'~','X');
        assignStructFieldsPerSpecies('n_sp',spnm{i},n,'cm^{-3}','n');
      end
      result = rmfield(result,'spgam');
    end

    % Create excited state number density fields
    if isfield(result,'exgam') && isfield(result,'rho')
      exspnm = fieldnames(result.exgam);
      for i = 1:numel(exspnm)
        % Array containing all of the excited state number densities
        ids = 1:min(length(result.exgam.(exspnm{i})),length(result.rho));
        n = result.exgam.(exspnm{i})(ids,:) .* result.rho(ids) .* avnCGS;
        assignStructFieldsPerSpecies('n_ex',exspnm{i},n,'cm^{-3}','n');
        % Assign subfields for each excited state
        fname = sprintf('n_ex_%s',exspnm{i});
        isp = find(strcmp(spnm,exspnm{i}));
        for iex = 1:mex(isp)
          en = levels.(exspnm{i}).energy(iex);
          [~,~,stateName] = getConfigTermField(exspnm{i},en);
          assignStructFieldsPerSpecies(fname,stateName,n(:,iex),'cm^{-3}','n');
        end
      end
      result = rmfield(result,'exgam');
    end
  end

  % Functions
  function assignStructFields(fieldName,resultValue,resultUnit,resultLabel)
    result.(fieldName) = resultValue;
    units.result.(fieldName)  = resultUnit;
    labels.result.(fieldName) = resultLabel;
  end
  function assignStructFieldsPerSpecies(fieldName,name,resultValue,resultUnit,resultLabel)
    result.(fieldName).(name) = resultValue;
    units.result.(fieldName).(name)  = resultUnit;
    labels.result.(fieldName).(name) = sprintf('%s_{%s}',resultLabel,name);
  end
  function assignStructFieldsAllSpecies(fieldName,resultValue,resultUnit,resultLabel)
    for s = 1:NS
      result.(fieldName).(spnm{s}) = resultValue(:,s);
      units.result.(fieldName).(spnm{s})  = resultUnit;
      if (resultLabel(end) == '}')
        labels.result.(fieldName).(spnm{s}) = sprintf('%s,%s}',resultLabel(1:end-1),spnm{s});
      else
        labels.result.(fieldName).(spnm{s}) = sprintf('%s_{%s}',resultLabel,spnm{s});
      end
    end
  end
  function data = fileData(latexVar)
    data = rawData(:,strcmp(latexVarNames,latexVar));
    if isempty(data)
      error('Variable %s not found in rawData by fileData.',latexVar);
    end
  end
end
