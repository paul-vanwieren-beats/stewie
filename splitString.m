function [fields] = splitString(line)

%Split line string into fields separated by commas
commas = [];
for i = 1:length(line)
    if line(i) == ','
        commas(end + 1) = i;
    end
end
if ~isempty(commas)
    n = length(commas);
    fields = cell(1, n + 1);
    fields{1} = line(1:(commas(1) - 1));
    for i = 2:n
        fields{i} = line((commas(i - 1) + 1):(commas(i) - 1));
    end
    fields{n + 1} = line((commas(end) + 1):end);
else
    fields = {line};
end