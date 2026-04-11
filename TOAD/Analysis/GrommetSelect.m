function Params = GrommetSelect(Index)
    % Read the text file as a comma-separated table
    T = readtable('GrommetOptions.txt', 'Delimiter', ',');
    
    % Simplify the field names to Name, Durometer, K, and B
    T.Properties.VariableNames = {'Name', 'Durometer', 'K', 'C'};
    
    % Validate the requested index
    if Index < 1 || Index > height(T)
        error('Index out of bounds. Please provide an index between 1 and %d.', height(T));
    end
    
    % Extract the requested row and convert it directly to a struct
    Params = table2struct(T(Index, :));
end