function structData = trimStructTERRA(structData,debug)
% Remove data at the end of every field within a struct to make everything
% the same length. Careful!

  if (nargin == 1)
    debug = true;
  end

  if isempty(fieldnames(structData))
    structData = structData;
    return;
  end

  % Get the minimum length within the struct
  minlen = 1e9;
  structFields = fieldnames(structData);
  nFields = max(size(structFields));
  flen = zeros(1,nFields);
  for rf = 1:nFields
    flen(rf) = fieldDataLength(structData.(structFields{rf}));
    minlen = min(flen(rf),minlen);
  end
  if debug
    fprintf("\nTrimming results to make them the same length. \n")
    for rf = 1:nFields
      fieldName = structFields{rf};
      fprintf("%s initial, final: \t %d \t %d \n", fieldName, flen(rf), minlen);
    end
  end

  % Stop at two layers deep, most expected within a TERRA struct
  trimIds = 1:minlen;
  mainFields = fieldnames(structData);
  for f1 = 1:numel(mainFields)
    fm = mainFields{f1};
    if isstruct(structData.(fm))
      secondaryFields = fieldnames(structData.(fm));
      for f2 = 1:numel(secondaryFields)
        fs = secondaryFields{f2};
        if isstruct(structData.(fm).(fs))
          error('Three layer deep structure is not supported for trimming!');
        else
          structData.(fm).(fs) = structData.(fm).(fs)(trimIds,:);
        end
      end
    else
      structData.(fm) = structData.(fm)(trimIds,:);
    end
  end

  if minlen == 0
    error('Minlen is zero in trimStructTERRA. Something is probably wrong!')
  end

  function len = fieldDataLength(structField)
    hasData = 0;
    while ~hasData
      if (isstruct(structField))
        fNames = fieldnames(structField);
        % Expect all fields to be same depth for result structure
        structField = structField.(fNames{1});
      else
        hasData = 1;
        len = length(structField(:,1));
      end
    end
  end
end
