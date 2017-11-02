
cd  '/Volumes/HD/Users/Marco/Documents/Projects/Matlab/BellmanFord'

maxi = realmax;

C  =     [     0   maxi   maxi   maxi   maxi
               3   maxi      4      2   maxi
            maxi      4   maxi      5      6
               7      2      5   maxi      4
            maxi   maxi      6      4   maxi ];
               
n   = size(C); 
v_1 = [0; maxi; maxi; maxi; maxi];
v   = v_1;

for t=1:n
	for j=1:n
		v(j) = min(v_1' + C(j,:));
	end;
	v_1 = v;
end;

 v'
		