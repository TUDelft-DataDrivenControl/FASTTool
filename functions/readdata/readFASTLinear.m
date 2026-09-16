function [data] = readFASTLinear(fileID)
% readFASTLinear extracts data of the linear model from FAST linearisation
% output. The extraction is separated by this function from opening the
% file to separate agnostic fileIO from the interpretation of the content
% of the file. Nevertheless, this function uses the fileID and reads
% directly from the file.
%
%   Syntax
%     [data] = readFASTLinear(fileID)
%
%   Input arguments
%     fileID - File identifier of an opened file with output from a FAST
%     linearisation
%
%   Output
%     data - Various data from the linear model described by the content of
%     the FAST linearisation output file
%
%   Called by
%     FileIOClass.readFASTLinOutput

% This function is called from FileIOClass.m, rather than from
% WindTurbineClass.m, which would be the better place to separate file
% management/reading from domain content. However, reading and interpreting
% this file is easier done when it can extract parts of the file, using its
% fileID.
%
% ToDo: This function and its helper functions need further commenting. The
% person who wrote it is no longer available for the development of
% FASTTool.

    % generic header:
                    fgetl(fileID); % skip a blank line
    data.ver{1,1} = fgetl(fileID); % FAST version info 
    data.ver{2,1} = fgetl(fileID); % submodule version info
                    fgetl(fileID); % skip a blank line
    data.desc     = fgetl(fileID); % model description
                    fgetl(fileID); % skip a blank line
                    fgetl(fileID); % "simulation information" comment/header line
                  
    % parse the next few lines:
    d = cell(8,1);
    for i=1:length(d)
        line = fgetl(fileID);
        C = textscan( line, '%s', 'delimiter', ':' );
        d{i} = textscan( C{1}{2}, '%f', 'CollectOutput',1 );
    end
    
    data.t        = d{1}{1};
    data.RotSpeed = d{2}{1};
    data.Azimuth  = d{3}{1};
    data.n_x      = d{4}{1};
    data.n_xd     = d{5}{1};
    data.n_z      = d{6}{1};
    data.n_u      = d{7}{1};
    data.n_y      = d{8}{1};
            

    line = fgetl(fileID);
    C = textscan( line, '%s', 'delimiter', '?' );
    if strfind( C{1}{2}, 'Yes' )
        SetOfMatrices = 2;
    else
        SetOfMatrices = 1;
    end         
    
    
    fgetl(fileID); % skip a blank line
    % get operating points and row/column order
    if data.n_x > 0 
        [data.x_op,    data.x_desc, data.x_rotFrame] = readLinTable(fileID,data.n_x);
        [data.xdot_op, data.xdot_desc]               = readLinTable(fileID,data.n_x);
    end

    if data.n_xd > 0 
        [data.xd_op,   data.xd_desc]                 = readLinTable(fileID,data.n_xd);
    end
    if data.n_z > 0 
        [data.z_op,    data.z_desc]                  = readLinTable(fileID,data.n_z);
    end
    if data.n_u > 0 
        [data.u_op,    data.u_desc, data.u_rotFrame] = readLinTable(fileID,data.n_u);
    end
    if data.n_y > 0 
        [data.y_op,    data.y_desc, data.y_rotFrame] = readLinTable(fileID,data.n_y);
    end
    
    
    fgetl(fileID); % skip a blank line
    for i=1:SetOfMatrices
        % get linearized state matrices
        fgetl(fileID); % skip linearized state matrices or jacobian description line
        fgetl(fileID); % skip a blank line

        while true
            [M, name] = readMatrix(fileID);
            if ~ischar(name) 
                break;
            end
            data.(name) = M;
        end
    end
end 

function [op, desc, RF] = readLinTable(fid,n)

    desc = cell(n,1);
    op   = cell(n,1);
    RF   = false(n,1);

    fgetl(fid); % table title/comment
    fgetl(fid); % table header row 1
    fgetl(fid); % table header row 2
    
    for row=1:n
        
        line = fgetl(fid);
        [C,pos] = textscan( line, '%*f %f %s',1 );
        if strcmp(line(pos),',') %we've got an orientation line:
            [C,pos] = textscan( line, '%*f %f %*s %f %*s %f %s',1 );
            op{row} = [C{1:3}];
        else
            op{row} = C{1};
        end
        RF(row) = strcmp(C{end},'T'); 
        desc{row}=strtrim( line(pos+1:end) );        
    end

    fgetl(fid); % skip a blank line
end

function [Mat, name] = readMatrix(fid)

    line = fgetl(fid);
    if ischar(line) && ~isempty(line)
        C = textscan( line, '%s', 'delimiter', ':' );
        name = C{1}{1};

        C = textscan( C{1}{2}, '%f %*s %f' );
        m=C{1};
        n=C{2};

        Mat = cell2mat( textscan(fid, repmat('%f',1,n),m,'CollectOutput',1) );
        fgetl(fid); %read end-of-line character(s)
    else
        name = -1;
        Mat = [];
    end
end
