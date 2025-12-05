function Link = ValveLink(name, type, ID, Up, Down, Rho, Cv, Active)
    Link.Name = name; Link.Type = type; Link.ID = ID; Link.Up = Up;
    Link.Down = Down; Link.Rho = Rho; Link.L = []; Link.A = []; Link.Zeta = [];
    Link.Cv = Cv; Link.Active = Active;
end