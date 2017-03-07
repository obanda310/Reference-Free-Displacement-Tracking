function newStack = reslice(imageStack,imagePlaneStr)
if imagePlaneStr == 'xz'
    newStack = flip(permute(imageStack,[3 2 1]),1);
end
if imagePlaneStr == 'yz'
    newStack = flip(permute(imageStack,[3 1 2]),1);
end
ShowStack(newStack)
end