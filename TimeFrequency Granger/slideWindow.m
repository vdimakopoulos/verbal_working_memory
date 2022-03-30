function indx = slideWindow(start, stop, winSize, shift)
% Create sliding windowed indices 

   
    
    startIndx = start : shift : stop;

    startIndx(find(startIndx >= stop)) = [];


    if (winSize == 1)
        stopIndx = startIndx + winSize;
    else
        stopIndx = startIndx + winSize - 1;
    end


    stopIndx(find(stopIndx > stop)) = stop;


    if (stopIndx(end) < stop)
        startIndx = horzcat(startIndx, stopIndx(end)+1);
        stopIndx = horzcat(stopIndx, stop);
    end

    indx = vertcat(startIndx, stopIndx);
end

