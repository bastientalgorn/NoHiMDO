function s = cleanSpaces(s)

s = strtrim(s);
s = deblank(s);

i = (s==' ');
i = find(i(1:end-1) & i(2:end));
s(i) = [];
