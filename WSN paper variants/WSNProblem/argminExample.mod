param N = 5;
param b{i in 1..N, j in 1..N} = i*j; 
param  minVal = min {i in 1..N, j in 1..N} b[i,j];

var x; var y;
for {i in 1..N, j in 1..N}{
	if b[i,j] = minVal then{
		let y := j;
		let x := i;
		break
	}
}

fix x;fix y;