function [qtot, qe, qv, qr] = partitionFunctionFromTERRA(T, levels, unitSystem)
  if nargin < 3
    unitSystem = 'SI';
  end
  switch unitSystem
    case 'SI'
      eVConversion = eV2J;
      boltz = boltzSI;
    case 'CGS'
      eVConversion = eV2erg;
      boltz = boltzCGS;
    otherwise
      error('Unrecognized unitSystem in partitionFunctionFromTERRA')
  end
  if iscell(T)
    if numel(T) ~= 3
      error('A cell array for T must have three entries for Trot, Tvib, Teex');
    end
    Tr = T{1};
    Tv = T{2};
    Tex = T{3};
    if length(Tr) ~= length(Tv) || length(Tv) ~= length(Tex)
      error('Temperature arrays must have same length.');
    end
  else
    Tr = T;
    Tv = T;
    Tex = T;
  end
  factr = levels.factr;
  ge = levels.g;
  ee = levels.energy * eVConversion;
  nex = length(ge);
  for iex = 1:nex
    eei = ee(iex);
    if ~isfield(levels, 'envx') % Atom
      qe(iex) = ge(iex) * exp(-eei / (boltz * T));
      qv = {};
      qr = {};
    else % Molecule
      ev{iex} = levels.envx{iex} * eVConversion;
      nv(iex) = length(ev{iex});
      for ivx = 1:nv(iex)
        er{iex}{ivx} = levels.enrx{ivx,iex} * eVConversion;
        evi = ev{iex}(ivx);
        nr{iex}(ivx) = length(er{iex}{ivx});
        for irx = 1:nr{iex}(ivx)
          eri = er{iex}{ivx}(irx);
          gr = levels.grv{ivx,iex}(irx);
          qr{iex}{ivx}(irx) = factr ...
                            * ge(iex) ...
                            * gr ...
                            * exp(-eei / (boltz * Tex)) ...
                            * exp(-evi / (boltz * Tv)) ...
                            * exp(-eri / (boltz * Tr));
        end
        qv{iex}(ivx) = sum(qr{iex}{ivx});
      end
      qe(iex) = sum(qv{iex});
    end
  end
  qtot = sum(qe);
end
