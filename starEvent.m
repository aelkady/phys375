function[value, isTerminal, direction] = starEvent(~, s)
global delta_t;
thresh = 1e-3;
if((delta_t <= thresh || (s(3) > 1.989e33)))
    value = 0;
else
    value = 1;
end;
isTerminal = 1;
direction = 0;