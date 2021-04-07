h = [-3 9 4 -10 9]
f = [6 -5 3]
y = myconv1d(h, f)
disp('Convolved = ')
disp(y)

function y = myconv1d(h, f)

    klen = length(f)
    padding = klen - 1

    hh = [zeros(1, padding) h zeros(1, padding)]
    ff = flip(f)

    y = zeros(1, length(h) + padding)
    for i = 1:length(hh)-padding
        y(i) = dot(ff, hh(i:i+padding))
    end
end