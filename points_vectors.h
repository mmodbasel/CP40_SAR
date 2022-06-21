/* 
 * File:   points_vectors.h
 * Author: martin
 *
 * Created on June 19, 2014, 5:02 PM
 */

#ifndef POINTS_VECTORS_H
#define	POINTS_VECTORS_H

/**********************************************/
/*                                            */
/*    P O I N T S   A N  D   V E C T O R S    */
/*                                            */
/**********************************************/

class Point {
	public:
	double	x, y, z;
	
	Point	operator+(Point p) { 
		Point p0;
		p0.x = x + p.x;
		p0.y = y + p.y;
		p0.z = z + p.z;
		return p0;
	}
	
	Point	operator-(Point p) { 
		Point p0;
		p0.x = x - p.x;
		p0.y = y - p.y;
		p0.z = z - p.z;
		return p0;
	}
};

class Vector {
	public:
	double	x, y, z;
	
	Vector	operator+ (Vector v) { 
		Vector v0;
		v0.x = x + v.x;
		v0.y = y + v.y;
		v0.z = z + v.z;
		return v0;
	}
	
	Vector	operator- (Vector v) { 
		Vector v0;
		v0.x = x - v.x;
		v0.y = y - v.y;
		v0.z = z - v.z;
		return v0;
	}
};

Vector	point2vector (Point p) {
	Vector	v;
	v.x = p.x;
	v.y = p.y;
	v.z = p.z;
	return v;
}

double	giveVectorLengthSquared(Vector v1) { return (v1.x*v1.x + v1.y*v1.y + v1.z*v1.z); }

double	giveDistance(Point A, Point B) {
        A = A - B;
        return sqrt(A.x*A.x + A.y*A.y + A.z*A.z);
}

double	giveDistanceSquared(Point A, Point B) {
        A = A - B;
        return (A.x*A.x + A.y*A.y + A.z*A.z);
}

double giveAngle(Point A, Point B, Point C) {
	double	asq, bsq, csq, arg;
	asq = giveVectorLengthSquared(point2vector(B - C));
        bsq = giveVectorLengthSquared(point2vector(A - C));
	csq = giveVectorLengthSquared(point2vector(B - A));
        arg = (bsq - asq - csq) / (-2.0 * sqrt(asq * csq));
        if(arg < -1.0) arg = -1.0;
        else if(arg > 1.0) arg = 1.0;
        return acos(arg);
}

/***********************/
/*                     */
/*    P H Y S I C S    */
/*                     */
/***********************/

// H-bonding routine based on the parameter by Angelo Vedani, Yeti Force field, JACS 1999
double	giveHbondEnergy(Point X, Point H, Point Y, unsigned char Nx, unsigned char Ny, int lowCutOff = NO){
	double	C, D;
	double	cosine;
	double	radius;
	double	wellDepth;
	double	rOpt;
	double	angleXHY;
        Point   p;
	
        p = H - Y;
        radius = p.x * p.x  +  p.y * p.y  +  p.z * p.z;
        
        
	// Arbitrary 3 A limit for H-bonds to O & N, 4 A to S
	if((Ny < SULPHUR  &&  radius < 9.0)  ||  radius < 16.0){
		angleXHY = giveAngle(X, H, Y);
		if(angleXHY > HB_THRESHOLD_ANGLE) {
                        //make sqrt out of dist
                        radius = sqrt(radius);
                        
			     if(Nx == 8  &&  Ny == 8)  { wellDepth = -4.946; rOpt = 1.746; }
			else if(Nx == 8  &&  Ny == 7)  { wellDepth = -4.655; rOpt = 1.878; }
			else if(Nx == 8  &&  Ny == 16) { wellDepth = -1.746; rOpt = 2.535; }
			else if(Nx == 7  &&  Ny == 8)  { wellDepth = -4.073; rOpt = 1.877; }
			else if(Nx == 7  &&  Ny == 7)  { wellDepth = -3.491; rOpt = 2.003; }
			else if(Nx == 7  &&  Ny == 16) { wellDepth = -1.455; rOpt = 2.667; }
			else if(Nx == 16 &&  Ny == 8)  { wellDepth = -2.328; rOpt = 2.099; }
			else if(Nx == 16 &&  Ny == 7)  { wellDepth = -2.037; rOpt = 2.088; }
			else if(Nx == 16 &&  Ny == 16) { wellDepth = -1.164; rOpt = 3.009; }
			else { printf("Error: Wrong atom type(s) (Nx = %d, Ny = %d)for hydrogen bonding passed.\n", Nx, Ny); exit(1); }

                        if(lowCutOff == YES  &&  radius < rOpt) radius = rOpt;

			C = -5 * wellDepth * pow(rOpt, 12);
			D = -6 * wellDepth * pow(rOpt, 10);

			cosine = cos(angleXHY);
			cosine = cosine * cosine;
			
			return ( cosine * (C/pow(radius, 12) - D/pow(radius, 10)) );
		}
	}
	return 0;
}

#endif	/* POINTS_VECTORS_H */

