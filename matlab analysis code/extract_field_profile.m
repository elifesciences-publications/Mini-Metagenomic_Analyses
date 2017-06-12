function output = extract_field_profile( structure_array, field_name, total_seq_length )
% Extracts the profile from a particular field of the entire array
% structure_array should be an array of structs
% field_name should be a char
% total_seq_length is a flag
% 2015.11.09 Brian Yu
% 2015.11.10 added functionality to extract profile of numeric fields

if isfield(structure_array, field_name)
    
    field_entry = cell(length(structure_array),1);
    
    if total_seq_length
        seqlength = zeros(size(structure_array));
    end
    
    for i = 1:length(structure_array)
        
        if ischar(structure_array(i).(field_name))
            field_entry{i} = structure_array(i).(field_name);
        elseif isnumeric(structure_array(i).(field_name)) && ...
                (length(structure_array(i).(field_name)) == 1)
            field_entry{i} = num2str(structure_array(i).(field_name));
        else
            error('Unknown type stored in entry %d field %s',i,field_name);
        end
        
        if total_seq_length && isfield(structure_array(i), 'contigLength')
            seqlength(i) = structure_array(i).contigLength;
        elseif ~isfield(structure_array(i), 'contigLength')
            error('Field contigLength does not exist in structure array');
        end
        
    end
    
    [C, ~, ic] = unique(field_entry);
    output = cell(length(C),2);
    if ~total_seq_length
        count = hist(ic,length(C));
        for i = 1:length(C)
            output{i,1} = count(i);
            output{i,2} = C{i};
        end
    else
        for i = 1:length(C)
            output{i,1} = 0;
            output{i,2} = C{i};
        end
        for i = 1:length(ic)
            output{ic(i),1} = output{ic(i),1} + structure_array(i).contigLength;
        end
    end
else
    output = nan;
end

end

