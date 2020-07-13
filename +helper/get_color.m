% darkgreen and purple go together
function out_color= get_color(col_name)

if ~exist('col_name', 'var')
    co= get(gca, 'colororder');
    col_num= get(gca, 'ColorOrderIndex');
    out_color= co(col_num+1, :);

elseif ischar(col_name)
    switch col_name
        case {'b', 'blue'}
            out_color= [69 117 180]/255;
            
        case {'r', 'red'}
            out_color= [215 48 39]/255;
            
        case {'g', 'green'}
            out_color= [27 158 119]/255;
            
        case {'lb', 'lightblue'}
            out_color= [145 191 219]/255;
            
        case {'lr', 'lightred'}
            out_color= [252 141 89]/255;
            
        case {'lg', 'lightgreen'}
            out_color= [102 194 165]/255;
            
        case {'dg', 'darkgreen'}
            out_color= [0 136 55]/255;

        case {'prp', 'purple'}
            out_color= [123 50 148]/255;
            
        case {'light_dg_prp'}
            out_color= [247 247 247]/255;
            
        case {'gray'}
            out_color= [100 100 100]/255;
            
        case {'black', 'k'}
            out_color= [0 0 0]/255;
        
        case {'brown', 'br'}
            out_color= [162 20 47]/255;
        
        case {'w', 'white'}
            out_color= [255 255 255]/255;
    end
    
elseif isnumeric(col_name)
    co= get(gca, 'colororder');
    out_color= co(col_name+1, :);
    
end