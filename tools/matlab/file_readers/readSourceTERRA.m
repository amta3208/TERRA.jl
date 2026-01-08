function [source, labels, units] = readSourceTERRA(units,labels,sourcePath,species,levels,prob_setup)
% This function reads source term output files from a TERRA run

  % Check if source terms were printed
  if isfield(prob_setup,'PRINT_SOURCE_TERMS') && (prob_setup.PRINT_SOURCE_TERMS == 0)
    source = struct();
    return
  end

  % Extract useful variables
  spnm = species.spnm;
  mex = species.mex;
  NS = length(mex);
  ies = species.ies;
  ih = species.ih;

  % Path stuff
  runPath = pwd;

  % Recognized source file names
  sourceNames = {'eeex','erot','evib','qrad','time'};
  if prob_setup.ND == 1
    sourceNames{end+1} = 'dist';
  end
  if isfield(prob_setup, 'PRINT_DISSOCIATION_SOURCES') && prob_setup.PRINT_DISSOCIATION_SOURCES == 1
    sourceNames{end+1} = 'diss';
  end
  if any(ies)
    sourceNames{end+1} = 'eion';
  end

  % Get the filenames
  sourceDir = dir(sourcePath);
  for i = 3:numel(sourceDir)
    sourceName = strrep(sourceDir(i).name,'.dat','');
    if strcmp2(sourceName,'species','start')
      sourceNames{end+1} = sourceName;
    end
  end

  % Read the source files and assign the data to the sources struct
  source = struct();
  for i = 1:numel(sourceNames)
    sourceType = sourceNames{i};
    sourceFile = sprintf('%s.dat',sourceType);
    sourceFilePath = sprintf('%s/%s',sourcePath,sourceFile);
    if ~isfile(sourceFilePath)
      error('Unrecognized file: %s', sourceFilePath);
    end
    rawData = importdata(sourceFilePath);
    switch sourceType
      case 'dist'
        [rawData, latexVarNames] = readOutputFileTERRA(sourcePath,sourceFile);
        assignStructFields('x',fileData('x [cm]'),'cm','x');
      case 'time'
        [rawData, latexVarNames] = readOutputFileTERRA(sourcePath,sourceFile);
        assignStructFields('t',fileData('t [\mus]'),'mu-s','t');
      case 'eeex'
        [rawData, latexVarNames] = readOutputFileTERRA(sourcePath,sourceFile);
        assignStructSubFields('eeex','eT',fileData('eT [erg/g-s]'),'erg/(g-s)','eT');
        assignStructSubFields('eeex','eR',fileData('eR [erg/g-s]'),'erg/(g-s)','eR');
        assignStructSubFields('eeex','eV',fileData('eV [erg/g-s]'),'erg/(g-s)','eV');
        assignStructSubFields('eeex','ECIE',fileData('E-CIE [erg/g-s]'),'erg/(g-s)','E-CIE');
        assignStructSubFields('eeex','ECII',fileData('E-CII [erg/g-s]'),'erg/(g-s)','E-CII');
        assignChemSources('eeex',rawData(:,find(strcmp(latexVarNames,'DISSOC [erg/g-s]')):end),'erg/(g-s)','e_{eex}');
      case 'erot'
        [rawData, latexVarNames] = readOutputFileTERRA(sourcePath,sourceFile);
        assignStructSubFields('erot','RT',fileData('RT [erg/g-s]'),'erg/(g-s)','RT');
        assignStructSubFields('erot','RV',fileData('RV [erg/g-s]'),'erg/(g-s)','RV');
        assignStructSubFields('erot','Re',fileData('Re [erg/g-s]'),'erg/(g-s)','Re');
        assignChemSources('erot',rawData(:,find(strcmp(latexVarNames,'DISSOC [erg/g-s]')):end),'erg/(g-s)','e_{rot}');
      case 'evib'
        [rawData, latexVarNames] = readOutputFileTERRA(sourcePath,sourceFile);
        assignStructSubFields('evib','VT',fileData('VT [erg/g-s]'),'erg/(g-s)','VT');
        assignStructSubFields('evib','VR',fileData('VR [erg/g-s]'),'erg/(g-s)','VR');
        assignStructSubFields('evib','Ve',fileData('Ve [erg/g-s]'),'erg/(g-s)','Ve');
        assignChemSources('evib',rawData(:,find(strcmp(latexVarNames,'DISSOC [erg/g-s]')):end),'erg/(g-s)','e_{vib}');
      case 'diss'
        [rawData, latexVarNames] = readOutputFileTERRA(sourcePath,sourceFile);
        for isp = 1:NS
          ispnm_true = strrep(strrep(spnm{isp},'p','+'),'n','-');
          if ~ies(isp) || (ih(isp) ~= 2)
            continue
          end
          for iex = 1:mex(isp)
            for jsp = 1:NS
              jspnm_true = strrep(strrep(spnm{jsp},'p','+'),'n','-');
              if strcmp(jspnm_true,'E')
                jspnm_true = 'E-';
              end
              fnm_eps_vib = sprintf('(\\langle\\epsilon_{vib}\\rangle_{d}/D_{0})_{%s(%d)-%s} [~]',ispnm_true,iex,jspnm_true);
              fnm_eps_ex  = sprintf('(\\langle\\epsilon_{ex}\\rangle_{d}/D_{0})_{%s(%d)-%s} [~]',ispnm_true,iex,jspnm_true);
              fnm_z = sprintf('Z_{%s(%d)-%s} [~]',ispnm_true,iex,jspnm_true);
              fnm_diss_rate = sprintf('fwd rate_{%s(%d)-%s} [mol/g-s]',ispnm_true,iex,jspnm_true);
              fnm_recomb_rate = sprintf('bwd rate_{%s(%d)-%s} [mol/g-s]',ispnm_true,iex,jspnm_true);
              fnm_net_diss_rate = sprintf('net rate_{%s(%d)-%s} [mol/g-s]',ispnm_true,iex,jspnm_true);
              assignStructSubFields(sprintf('eps_vib_diss_%s_%d',spnm{isp},iex),spnm{jsp},fileData(fnm_eps_vib),'~',strrep(fnm_eps_vib,' [~]',''))
              assignStructSubFields(sprintf('eps_ex_diss_%s_%d',spnm{isp},iex),spnm{jsp},fileData(fnm_eps_ex),'~',strrep(fnm_eps_ex,' [~]',''))
              assignStructSubFields(sprintf('Z_%s_%d',spnm{isp},iex),spnm{jsp},fileData(fnm_z),'~',strrep(fnm_z,' [~]',''))
              assignStructSubFields(sprintf('diss_rate_%s_%d',spnm{isp},iex),spnm{jsp},fileData(fnm_diss_rate),'~',strrep(fnm_diss_rate,' [mol/g-s]',''))
              assignStructSubFields(sprintf('recomb_rate_%s_%d',spnm{isp},iex),spnm{jsp},fileData(fnm_recomb_rate),'~',strrep(fnm_recomb_rate,' [mol/g-s]',''))
              assignStructSubFields(sprintf('net_rate_%s_%d',spnm{isp},iex),spnm{jsp},fileData(fnm_net_diss_rate),'~',strrep(fnm_net_diss_rate,' [mol/g-s]',''))
            end
          end
        end
      case 'eion'
        [rawData, latexVarNames] = readOutputFileTERRA(sourcePath,sourceFile);
        for isp = 1:NS
          if ~ies(isp) || (ih(isp) == 0)
            continue
          end
          ispnm_true = strrep(strrep(spnm{isp},'p','+'),'n','-');
          fnm_alpha = sprintf('(\\langle\\epsilon_{eex}\\rangle_{d}/I_{0})_{%s} [~]',ispnm_true);
          fnm_z_eion = sprintf('Z_{%s, EII} [~]',ispnm_true);
          assignStructSubFields('alpha_eii',spnm{isp},fileData(fnm_alpha),'~',strrep(fnm_alpha,' [~]',''))
          assignStructSubFields('Z_eii',spnm{isp},fileData(fnm_z_eion),'~',strrep(fnm_z_eion,' [~]',''))
        end
      case 'qrad'
        [rawData, latexVarNames] = readOutputFileTERRA(sourcePath,sourceFile);
        assignStructSubFields('qrad','radbb',fileData('bound-bound [erg/g-s]'),'erg/(g-s)','bound-bound');
        assignStructSubFields('qrad','radfb',fileData('free-bound [erg/g-s]'),'erg/(g-s)','free-bound');
        assignStructSubFields('qrad','radff',fileData('free-free [erg/g-s]'),'erg/(g-s)','free-free');
      otherwise
        if contains(sourceType,'species')
          [rawData, latexVarNames] = readOutputFileTERRA(sourcePath,sourceFile);
          if contains(sourceType,'state')
            ids = sscanf(sourceType,'species_%d_state_%d');
            spid = ids(1);
            exid = ids(2);
            speciesName = spnm{spid};
            levelEnergy = levels.(speciesName).energy(exid);
            [~,stateTerm,stateName] = getConfigTermField(speciesName,levelEnergy);
            nameWithTerm = sprintf('%s(%s)',speciesName,stateTerm);
            assignChemSources(stateName,rawData(:,1:find(strcmp(latexVarNames,'Other  [mol/g-s]'))),'mol/(g-s)',nameWithTerm);
            assignStructSubFields(stateName,'HCIE',fileData('H-CIE [mol/g-s]'),'mol/(g-s)',sprintf('%s, H-CIE',nameWithTerm));
            assignStructSubFields(stateName,'HCII',fileData('H-CII [mol/g-s]'),'mol/(g-s)',sprintf('%s, H-CII',nameWithTerm));
            assignStructSubFields(stateName,'ECIE',fileData('E-CIE [mol/g-s]'),'mol/(g-s)',sprintf('%s, E-CIE',nameWithTerm));
            assignStructSubFields(stateName,'ECII',fileData('E-CII [mol/g-s]'),'mol/(g-s)',sprintf('%s, E-CII',nameWithTerm));
            assignStructSubFields(stateName,'radbb',fileData('bound-bound radiation [mol/g-s]'),'mol/(g-s)',sprintf('%s, bound-bound radiative transition',nameWithTerm));
            assignStructSubFields(stateName,'radfb',fileData('free-bound radiation [mol/g-s]'),'mol/(g-s)',sprintf('%s, radiative/dielectronic recombination',nameWithTerm));
          else
            spid = sscanf(sourceType,'species_%d');
            speciesName = spnm{spid};
            assignChemSources(speciesName,rawData(:,1:find(strcmp(latexVarNames,'Other  [mol/g-s]'))),'mol/(g-s)',speciesName);
            assignStructSubFields(speciesName,'HCII',fileData('H-CII [mol/g-s]'),'mol/(g-s)',sprintf('%s, H-CII',speciesName));
            assignStructSubFields(speciesName,'ECII',fileData('E-CII [mol/g-s]'),'mol/(g-s)',sprintf('%s, E-CII',speciesName));
            assignStructSubFields(speciesName,'radfb',fileData('free-bound radiation [mol/g-s]'),'mol/(g-s)',sprintf('%s, radiative/dielectronic recombination',speciesName));
          end
        else
          error('Unrecognized source file type');
        end
    end
  end

  cd(runPath);

  % Functions
  function assignStructFields(fieldName,sourceValue,sourceUnit,sourceLabel)
    source.(fieldName) = sourceValue;
    units.source.(fieldName)  = sourceUnit;
    labels.source.(fieldName) = sourceLabel;
  end
  function assignStructSubFields(fieldName,subFieldName,sourceValue,sourceUnit,sourceLabel)
    source.(fieldName).(subFieldName) = sourceValue;
    units.source.(fieldName).(subFieldName)  = sourceUnit;
    labels.source.(fieldName).(subFieldName) = sourceLabel;
  end
  function assignChemSources(fieldName,sourceValues,sourceUnits,sourceFieldName)
    if length(sourceValues(1,:)) ~= 7
      error('Received more than 7 variables to assign. Check source file formatting against the expected formatting implemented here.')
    end
    assignStructSubFields(fieldName,'DISSOC',sourceValues(:,1),sourceUnits,sprintf('%s, Dissociation',sourceFieldName));
    assignStructSubFields(fieldName,'NEUTEX',sourceValues(:,2),sourceUnits,sprintf('%s, Neutral Exchange',sourceFieldName));
    assignStructSubFields(fieldName,'AI_DR',sourceValues(:,3),sourceUnits,sprintf('%s, AI/DR',sourceFieldName));
    assignStructSubFields(fieldName,'CHAREX',sourceValues(:,4),sourceUnits,sprintf('%s, Charge Exchange',sourceFieldName));
    assignStructSubFields(fieldName,'ELEXNG',sourceValues(:,5),sourceUnits,sprintf('%s, Excitation Exchange',sourceFieldName));
    assignStructSubFields(fieldName,'QUENCH',sourceValues(:,6),sourceUnits,sprintf('%s, Quenching',sourceFieldName));
    assignStructSubFields(fieldName,'OtherChem',sourceValues(:,7),sourceUnits,sprintf('%s, Other Chemistry',sourceFieldName));
  end
  function data = fileData(latexVar)
    data = rawData(:,strcmp(latexVarNames,latexVar));
    if isempty(data)
      error('Variable %s not found in rawData by fileData.',latexVar);
    end
  end
end
