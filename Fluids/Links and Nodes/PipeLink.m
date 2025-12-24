function Link = PipeLink(name, ID, Up, Down, L, A, Zeta)
Link.Name = name; Link.Type = 'Pipe'; Link.ID = ID; Link.Up = Up;
Link.Down = Down; Link.L = L; Link.A = A; Link.Zeta = Zeta; 
Link.Cv = [];
end