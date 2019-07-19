%Create data from a function corrupted with gaussian noise
function data = generate_data(func,noise,pts)
    data = zeros(pts,2);
    for i = 1:pts
        data(i,1) = rand();
        data(i,2) = feval(func,data(i,1)) + sqrt(noise)*randn();
    end
end