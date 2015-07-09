function [decimal] = hexadecimalToDecimal(string)

%Convert hexadecimal number to decimal number
n = length(string);
idx = 0;
decimal = 0;
for i = n:(-1):1
    
    switch string(i)
        case '1'
            decimal = decimal + 1*16^idx;
        case '2'
            decimal = decimal + 2*16^idx;
        case '3'
            decimal = decimal + 3*16^idx;
        case '4'
            decimal = decimal + 4*16^idx; 
        case '5'
            decimal = decimal + 5*16^idx; 
        case '6'
            decimal = decimal + 6*16^idx; 
        case '7'
            decimal = decimal + 7*16^idx; 
        case '8'
            decimal = decimal + 8*16^idx; 
        case '9'
            decimal = decimal + 9*16^idx; 
        case {'A', 'a'}
            decimal = decimal + 10*16^idx;
        case {'B', 'b'}
            decimal = decimal + 11*16^idx; 
        case {'C', 'c'}
            decimal = decimal + 12*16^idx; 
        case {'D', 'd'}
            decimal = decimal + 13*16^idx;
        case {'E', 'e'}
            decimal = decimal + 14*16^idx; 
        case {'F', 'f'}
            decimal = decimal + 15*16^idx; 
    end
    
    idx = idx + 1;
    
end

%Convert to signed
iMax = 2^(idx * 4);
iNeg = (iMax / 2) - 1;
decimal = decimal - iMax*(decimal > iNeg);
    