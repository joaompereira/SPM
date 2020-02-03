function value = setordefault(options, field, default_value)
% Merge input options with local defaults

  if isfield(options,field)
      value = options.(field);
  else
      value = default_value;
  end
  
end