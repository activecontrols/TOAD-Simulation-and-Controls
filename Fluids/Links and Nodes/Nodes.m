function Node = Nodes(name, ID, Fixed, P0, V, Y0, LinksIN, LinksOUT, isCombustor)

Node(1).Name = name; Node(1).ID = ID; Node(1).Fixed = Fixed; Node(1).P0 = P0;
Node.V = V; Node.Y0 = Y0; Node(1).LinksIN = LinksIN; Node(1).LinksOUT = LinksOUT;
if isCombustor
    CstarEff = 0.83;
    Node.Temp = 2666 * CstarEff^2; Node.R = 442; Node.Gamma = 1.207; Node.IsCombustor = true;
else
    Node.Temp = 293; Node.R = 296.8; Node.Gamma = 1.2; Node.IsCombustor = false;
end