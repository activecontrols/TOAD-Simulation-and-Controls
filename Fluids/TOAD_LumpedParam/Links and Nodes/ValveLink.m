function Link = ValveLink(name, type, ID, Up, Down, Cv, varargin)
    Link.Name = name; Link.Type = type; Link.ID = ID; Link.Up = Up;
    Link.Down = Down; Link.L = []; Link.A = 0; Link.Zeta = [];
    Link.Cv = Cv;
    
    % Pre-allocate Regulator fields for struct uniformity
    Link.P_set = 0;
    Link.Droop = 1; 
    Link.SPE   = 0;
    
    if strcmp(type, 'Orifice')
        if ~isempty(varargin)
            Link.A = varargin{1};
        else
            error('Error: Link "%s" is an Orifice but no Area was provided!', name);
        end
    elseif strcmp(type, 'Regulator')
        if length(varargin) >= 3
            Link.P_set = varargin{1} * 6895; 
            Link.Droop = varargin{2} * 6895; 
            Link.SPE   = varargin{3};
        else
            error('Error: Link "%s" is a Regulator. Provide P_set, Droop, and SPE.', name);
        end
    end
end