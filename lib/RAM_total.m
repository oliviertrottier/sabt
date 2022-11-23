function Mem = RAM_total()
% Function to calculate the available memory (GB) for Unix and Windows
% OS.
parse_num = @(s) str2double(regexp(s,'\d+','match','once'));
if ismac()
    % Use vm_stat to query available RAM.
    [~,Res] = system('vm_stat');
    Res_lines = strsplit(Res,'\n');
    
    % Determine the page size in bytes.
    if contains(Res_lines{1},'page size')
        Page_size = regexp(Res_lines{1},'\d+','match');
        if numel(Page_size) > 1
            warning('Cannot infer page size. Using 4096 bytes.');
            Page_size = 4096;
        else
            Page_size = str2double(Page_size{1});
        end
    else
        Page_size = 4096;
    end
    
    % Parse each line to calculate the total memory.
    Mem = 0;
    for i=2:numel(Res_lines)
        if contains(Res_lines{i},{'Pages free','Pages active','Pages inactive','Pages wired down'})
            Mem = Mem + parse_num(Res_lines{i});
        end
    end
    Mem = Mem * Page_size/(1024^3); % Memory in GB

elseif isunix()
    % Use vmstat to query available RAM.
    [~,Res] = system('vmstat --stats --unit K'); % Query memory in kB.
    Res_lines = strsplit(Res,'\n');
    
    % Find the line with the total memory
    Mem = 0;
    for i=1:numel(Res_lines)
        if contains(Res_lines{i},{'total memory'})
            Mem = Mem + parse_num(Res_lines{i});
        end
    end
    Mem = Mem/(1024^2); % Memory in GB
else
    % Windows. Use the built-in memory command.
    [~,systemview] = memory();
    Mem = systemview.PhysicalMemory.Total/(1024^3);
end

% Check that the memory was properly calculated.
assert(Mem > 0,'The total RAM could not be calculated.');
end