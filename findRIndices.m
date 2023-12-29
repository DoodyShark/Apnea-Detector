% This function finds the indeces at which the R peaks were found
function RIndeces = findRIndices(signal)
    RIndeces = [];
    for i = 2:length(signal)-1
        if ((signal(i) > signal(i-1)) && (signal(i) >= signal(i+1)) && (signal(i) > 0.6))
            RIndeces = [RIndeces, i];
        end
    end
end