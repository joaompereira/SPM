function varargout = option_parser(user_options, varargin)
% Merge input options with local defaults
  
  n = length(varargin);
  iP = inputParser;
  % vargin are the default options
  for i=1:n
      addParameter(iP,varargin{i}{:});
  end
  parse(iP,user_options{:});
  
  varargout = cell(1,n);
  for i=1:length(varargin)
      varargout{i} = iP.Results.(varargin{i}{1});
  end
  
end