function poo = testfunhandlesmaster()
%     a = 2;
    fun = @helperfun;
    
    for j = 1:5
        a = j;
        poo(j) = testfunhandlesaux(fun);
    end
    
    function O = helperfun(x)
        O = a*x^2;
    end

end