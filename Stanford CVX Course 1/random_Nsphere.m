% Returns a random vector lying on surface 'n' dimensional unit sphere
function x = random_Nsphere(n)
    x = rand(1,n);
    
    if n<1
        fprintf('N should be >= 1');
        x = NaN;
    else
        if n > 1          
            sum = 0;
            for i = 1:n
               sum = sum + x(1,i)*x(1,i);
            end
            sum = sum^0.5;
            x = x./sum;
        end
    end
end