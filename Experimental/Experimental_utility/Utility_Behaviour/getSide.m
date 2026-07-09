%function side = getSide(bias)

function side = getSide(bias)
if bias > 0
    side = "CW";
else 
    side = "CCW";
end
end