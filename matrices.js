/*
   matrices.js:
   function to handle matrices
   Courtesy of Marco Tarini, Università dell'Insubria
*/

/* ALL MATRICES ARE IN COLUMN-MAJOR FORMAT  */
/* ALL ANGLES ARE IN DEGREES                */

function identityMatrix() {
	return [
		1,0,0,0,
		0,1,0,0,
		0,0,1,0,
		0,0,0,1
	];
}

function rotationXMatrix( deg ){
	var s = Math.sin( deg * 3.1415 / 180);
	var c = Math.cos( deg * 3.1415 / 180);
	return res = [
	   1 ,0 ,0 ,0,
	   0 ,c ,s ,0,
	   0 ,-s,c ,0,
	   0 ,0, 0 ,1
	];
}

function rotationYMatrix( deg ){
	var s = Math.sin( deg * 3.1415 / 180);
	var c = Math.cos( deg * 3.1415 / 180);
	return res = [
	   c ,0 ,-s,0,
	   0 ,1 ,0 ,0,
	   s ,0, c ,0,
	   0 ,0, 0 ,1
	];
}

function rotationZMatrix( deg ){
	var s = Math.sin( deg * 3.1415 / 180);
	var c = Math.cos( deg * 3.1415 / 180);
	return res = [
	   c ,s ,0 ,0,
	   -s,c ,0 ,0,
	   0 ,0, 1 ,0,
	   0 ,0, 0 ,1
	];
}

function translationMatrix( dx, dy, dz ){
	return [
	   1 ,0 ,0 ,0,
	   0 ,1 ,0 ,0,
	   0 ,0 ,1 ,0,
	   dx,dy,dz,1
	];
}

function scalingMatrix( sx, sy, sz ){
    if (typeof sy == 'undefined' ) sy = sx;
    if (typeof sz == 'undefined' ) sz = sx;
	return [
	   sx,0 ,0 ,0,
	   0 ,sy,0 ,0,
	   0 ,0 ,sz,0,
	   0 ,0 ,0 ,1
	];
}

function perspectiveMatrix( focal, aspect, near, far ){
	return [
	   focal ,0           ,0 ,                     0,
	   0     ,focal*aspect,0 ,                     0,
	   0     ,0           ,(far+near)/(near-far) ,-1,
	   0     ,0           ,2*far*near/(near-far) , 0
	];
}

function perspectiveMatrixFOV( fov, aspect, near, far ){
    var focal = 1.0/Math.tan(fov * 3.1415 / 360);
	return [
	   focal/aspect,   0,     0,                      0,
	   0           ,   focal, 0,                      0,
	   0           ,   0,     (far+near)/(near-far), -1,
	   0           ,   0,     2*far*near/(near-far),  0
	];
}

function multMatrix( a, b ) {
	/* product row by column */
	var res = [];
	for (var i=0; i<4; i++)    // column
	for (var j=0; j<4; j++) {  // row
		res[i*4+j] = 0;
		for (var k=0; k<4; k++)
			res[i*4+j] += a[k*4+j] * b[i*4+k];
	}
	return res;
}

function multMatrixVec( a, v ) {
    var res = [];
    for (var j=0; j<4; j++) {  // row
        res[j]=0;
        for (var k=0;k<4;k++)
            res[j] += a[k*4+j]*v[k];
    }
    return res;
}

function normVec( v ) {
    var n = Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    return [v[0]/n,v[1]/n,v[2]/n,0]
}

function transposeMatrix3(m) {
    return [m[0], m[3], m[6],
            m[1], m[4], m[7],
            m[2], m[5], m[7]];
}

function transposeMatrix4(m) {
    return [m[0], m[4], m[ 8], m[12],
            m[1], m[5], m[ 9], m[13],
            m[2], m[6], m[10], m[14],
            m[3], m[7], m[11], m[15]];
}


function detMatrix3(m) {
    return  m[0]*(m[4]*m[8]-m[5]*m[7]) -
            m[1]*(m[3]*m[8]-m[5]*m[6]) +
            m[2]*(m[3]*m[7]-m[4]*m[6]);
}

function invMatrix3(m) {
    var det = detMatrix3(m);
    if (det == 0) return [0,0,0,0,0,0,0,0,0];
    var A = m[4]*m[8]-m[7]*m[5];
    var B = m[7]*m[2]-m[1]*m[8];
    var C = m[5]*m[1]-m[4]*m[2];
    var D = m[6]*m[5]-m[3]*m[8];
    var E = m[0]*m[8]-m[6]*m[2];
    var F = m[3]*m[2]-m[0]*m[5];
    var G = m[3]*m[7]-m[6]*m[4];
    var H = m[6]*m[1]-m[0]*m[7];
    var I = m[4]*m[0]-m[3]*m[1];
    return [A/det, B/det, C/det,
            D/det, E/det, F/det,
            G/det, H/det, I/det];
}

function invMatrix4(m) {
    var a = [m[0], m[1], m[2],
             m[4], m[5], m[6],
             m[8], m[9], m[10]];
    var b = [m[12], m[13], m[14]];
    a = invMatrix3(a);
    var c = [-a[0]*b[0]-a[3]*b[1]-a[6]*b[2],
             -a[1]*b[0]-a[4]*b[1]-a[7]*b[2],
             -a[2]*b[0]-a[5]*b[1]-a[8]*b[2]];
    return [a[0], a[1], a[2], 0,
            a[3], a[4], a[5], 0,
            a[6], a[7], a[8], 0,
            c[0], c[1], c[2], 1];
}

/* Vector utilities */
function vec3( x,y,z ){
	return [x,y,z];
}

function subtract( a,b ){
	return [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ];
}

function sum( a,b ){
	return [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ];
}

function normalize( a ){
	k = Math.sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
	return [ a[0]/k, a[1]/k, a[2]/k ];
}

//additional functions
function multScalarVector(s, v) {
	var result = [];
	for(var i = 0; i < v.length; i++)
		result[i] = v[i] * s;
	return result;
}

function crossProduct(v0, v1) {
	return [v0[1]*v1[2] - v0[2]*v1[1],
          	v0[2]*v1[0] - v0[0]*v1[2],
          	v0[0]*v1[1] - v0[1]*v1[0]];
}

function magnitudeVector(v) {
	var magnitude = 0;
	for(var i = 0; i < v.length; i++)
		magnitude += Math.pow(v[i], 2);
	if(magnitude == 0)
		return magnitude;
	return Math.sqrt(magnitude);
}

function sumVectors(v1, v2, v3) {
	if(v3 == null)
		return [v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]]
	return [v1[0]+v2[0]+v3[0], v1[1]+v2[1]+v3[1], v1[2]+v2[2]+v3[2]]
}

function dotProduct(v1, v2) {
	var dotprod = 0;
	for(var i = 0; i < 3; i++)
		dotprod += v1[i] * v2[i];
	return dotprod;
}

//quaternion functions
function identityQuat() {
	return [[0, 0, 0], 1];
}

function axisAngle2Quat(axis, angle) {
	var	w = Math.cos(angle/2),
		sin = Math.sin(angle/2);

	return [[axis[0]*sin, axis[1]*sin, axis[2]*sin], w];
}

function multQuat(q1, q2) {
	return [sumVectors(multScalarVector(q2[1],q1[0]),
					  multScalarVector(q1[1], q2[0]),
					  crossProduct(q1[0], q2[0])),
		  	q1[1]*q2[1]-dotProduct(q1[0], q2[0])];
}

function quat2Matrix(q) {
	var x = q[0][0],
		y = q[0][1],
		z = q[0][2],
		w = q[1];

	//adding a column and a row because we always deal with 4x4 matrices
	return [
		1-2*Math.pow(y, 2)-2*Math.pow(z, 2), 2*x*y-2*w*z, 2*x*z+2*w*y, 0,
		2*x*y+2*w*z, 1-2*Math.pow(x, 2)-2*Math.pow(z, 2), 2*x*y-2*w*x, 0,
		2*x*z-2*w*y, 2*y*z+w*x, 1-2*Math.pow(x, 2)-2*Math.pow(y, 2), 0,
		0, 0, 0, 1
	]
}

function magnitudeQuat(q) {
	var x = q[0][0],
		y = q[0][1],
		z = q[0][2],
		w = q[1];

	return Math.sqrt(Math.pow(w, 2) + Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2))
}

function normQuat(q) {
	return q/magnitudeQuat(q);
}
