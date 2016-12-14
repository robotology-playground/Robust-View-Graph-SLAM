% c = setminus(a,b)  Set-theoretical subtraction of index sets a, b.

function c = setminus(a,b);
c = a(~ismembc(a,sort(double(b))));
return
