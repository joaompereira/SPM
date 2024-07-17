function opts = option_parser(user_options, varargin)
% Merge input options with local defaults
  
  n = length(varargin);
  iP = inputParser;
  % vargin are the default options
  for i=1:n
      addParameter(iP,varargin{i}{:});
  end
  parse(iP,user_options{:});
  
  for i=1:n
      opts.(varargin{i}{1}) = iP.Results.(varargin{i}{1});
  end
  
end