% repressive hill function,  50% level for repression is sat
function rep_hill = repression(species, sat)
rep_hill = sat/(sat+species);
end
