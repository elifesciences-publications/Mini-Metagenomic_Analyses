function output = extract_structarray_entries( structure_array, field_name, value )
% Extracts the profile from a particular field of the entire array
% structure_array should be an array of structs
% field_name should be a char
% value is the identity of the field entry ie. 'Bacteria'
% 2015.11.09 Brian Yu
% 2015.11.10 Added functionality so that entry field could be numbers

if isfield(structure_array, field_name)
    selid = false(size(structure_array));
    for i = 1:length(structure_array)
        if ischar(structure_array(i).(field_name)) && strcmpi(structure_array(i).(field_name), value)
            selid(i) = true;
        elseif isnumeric(structure_array(i).(field_name)) && ...
                (length(structure_array(i).(field_name)) == 1) && ...
                strcmpi(num2str(structure_array(i).(field_name)), num2str(value))
            selid(i) = true;
        end
    end
    output = structure_array(selid);
else
    output = nan;
end

end

