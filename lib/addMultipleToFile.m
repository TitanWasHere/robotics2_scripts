function addMultipleToFile(varargin)
% addMultipleToFile Appends multiple input parameters to output.txt.
%
%   addMultipleToFile(PARAM1, PARAM2, ...) opens (or creates) output.txt in 
%   append mode and writes each parameter's string representation on a new line.
%
% Example:
%   addMultipleToFile(123, 'hello', [1 2 3]);

    % Open file for appending
    fid = fopen('output.txt', 'a');
    if fid == -1
        error('Error opening output.txt for writing.');
    end

    % Loop over each input parameter
    for i = 1:nargin
        param = varargin{i};

        % Convert parameter to string based on its type
        if isnumeric(param)
            paramStr = num2str(param);
        elseif ischar(param) || isstring(param)
            paramStr = char(param); % Convert string to character vector if needed
        else
            % Use JSON encoding for other data types
            paramStr = jsonencode(param);
        end

        % Write the parameter to the file, followed by a newline.
        fprintf(fid, '%s\n', paramStr);
    end
    
    % Close the file.
    fclose(fid);
end
