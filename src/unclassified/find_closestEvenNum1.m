function num = find_closestEvenNum1(numIn)
    
    tn = round(numIn * 100);
    if mod(tn, 2) == 0, 
        tn = tn; 
    else
        tn = tn + 1;
    end
    
    num = tn / 100;

