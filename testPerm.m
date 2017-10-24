%%

numTest = 1000;
NP = 320;
rot = (0:1:NP-1);


for cnt = 1:numTest
  ind = randperm(NP-1);              % index pointer array
  
  a1  = randperm(NP);             % shuffle locations of vectors
  rt2 = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
  a2  = a1(rt2+1);                 % rotate vector locations
  rt3 = rem(rot+ind(2),NP);
  a3  = a1(rt3+1);
  rt4 = rem(rot+ind(3),NP);
  a4  = a1(rt4+1);
  rt5 = rem(rot+ind(4),NP);
  a5  = a1(rt5+1);
  
  failTest = any(a1==a2) | any(a1==a3) | any(a1==a4) | any(a1==a5) | any(a2==a3)...
  | any(a2==a4) | any(a2==a5) | any(a3==a4) | any(a3==a5) | any(a4==a5);
  if failTest
    error('You Failed The Test')
  end
endfor
