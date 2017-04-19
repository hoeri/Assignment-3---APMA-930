%% Test saving data in a script

function    SaveTest

N = 20 ;

Squares = zeros(1,N) ;

for n = 1:N
    Squares(n) = n.^2 ;
end

FileName = 'Squares.mat' ;

save(FileName,'Squares')  ;





