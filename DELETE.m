function [x] = test_func(x, p)
    p = inputParser;
    argName = 'myInput';
    defaultVal = 13;
    addOptional(p, argName, defaultVal)
end
