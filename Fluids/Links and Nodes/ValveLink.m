function Link = ValveLink(name, type, ID, Up, Down, Rho, Cv, varargin)
    Link.Name = name; Link.Type = type; Link.ID = ID; Link.Up = Up;
    Link.Down = Down; Link.Rho = Rho; Link.L = []; Link.A = 0; Link.Zeta = [];
    Link.Cv = Cv;
    if strcmp(type, 'Orifice')
        if ~isempty(varargin)
            Link.A = varargin{1};
        else
            error('Error: Link "%s" is an Orifice but no Area was provided!', name);
        end
    end
end