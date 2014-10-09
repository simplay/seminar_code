function t = intervalSample(n)


intervalB = 0:(1.0/(n+1)):1;

t = zeros(n,1);
for k=1:n,
    t(k) = rand(1,1)*(intervalB(k+1)-intervalB(k))+intervalB(k+1);
end

end