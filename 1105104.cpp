#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include"bitmap_image.hpp"
using namespace std;

#define pi (2*acos(0.0))
#define EPSILON 0.001
#define EPSILON_SQUARE (EPSILON*EPSILON)
#define FAR 1.0
#define CHECK <
#define zBuf +

//Vector and Point class definitions below
class Vector;
class Point
{
public:
	double arr[4];
	double getX()
	{
	    return arr[0];
	}
	double getY()
	{
	    return arr[1];
	}
	double getZ()
	{
	    return arr[2];
	}
	Point()
	{
	    for(int i = 0; i<3; i++)
        {
            arr[i] = 0;
        }
        arr[3] = 1;
	}
	Point(double x, double y, double z, double w)
	{
	    arr[0] = x/w;
	    arr[1] = y/w;
	    arr[2] = z/w;
	    arr[3] = 1;
	}
	Point& operator=(const Point& p)
	{
	    //cout << "assigning point" << endl;
	    this->arr[0] = p.arr[0];
	    this->arr[1] = p.arr[1];
	    this->arr[2] = p.arr[2];
	    this->arr[3] = p.arr[3];
	    return *this;
	}
    Point operator+(const Vector& v);
	Vector operator-(const Point& p);
	void print()
	{
	    cout << "(" << arr[0] << ", " << arr[1] << ", " << arr[2]  << ")" << endl;
	}
};

class Vector
{
public:
	double x, y, z;
	Vector()
	{
	    x = 0;
	    y = 0;
	    z = 0;
	}
	Vector(double x, double y, double z)
	{
	    this->x = x;
	    this->y = y;
	    this->z = z;
	}
    void print()
	{
	    cout << "[" << x << ", " << y << ", " << z << "]" << endl;
	}
	Vector operator+(const Vector& v)
	{
		Vector retV;
		retV.x = this->x + v.x;
		retV.y = this->y + v.y;
		retV.z = this->z + v.z;
		return retV;
	}
	Point operator+(const Point& p)
	{
		Point retV;
		retV.arr[0] = this->x + p.arr[0]/p.arr[3];
		retV.arr[1] = this->y + p.arr[1]/p.arr[3];
		retV.arr[2] = this->z + p.arr[2]/p.arr[3];
		retV.arr[3] = 1;
		return retV;
	}
	Vector operator-(const Vector& v)
	{
		Vector retV;
		retV.x = this->x - v.x;
		retV.y = this->y - v.y;
		retV.z = this->z - v.z;
		return retV;
	}
	Vector operator*(const Vector& v)
	{
		Vector retValue;
		retValue.x = this->y*v.z - this->z*v.y;
		retValue.y = this->z*v.x - this->x*v.z;
		retValue.z = this->x*v.y - this->y*v.x;
		return retValue;
	}
	double operator,(const Vector& v)
	{
		double retV;
		retV = this->x * v.x + this->y * v.y + this->z * v.z;
		return retV;
	}
	Vector operator*(double s)
	{
		Vector retV;
		retV.x = this->x * s;
		retV.y = this->y * s;
		retV.z = this->z * s;
		return retV;
	}
	void normalize()
	{
	    double magnitude = (this->x * this->x) + (this->y * this->y) + (this->z * this->z);
	    magnitude = sqrt(magnitude);
	    //cout << "magnitude = " << magnitude << endl;
	    this->x = this->x / magnitude;
	    this->y = this->y / magnitude;
	    this->z = this->z / magnitude;
	}
};

Point Point::operator+(const Vector& v)
{
	Point retV;
	retV.arr[0] = this->arr[0]/this->arr[3] + v.x;
	retV.arr[1] = this->arr[1]/this->arr[3] + v.y;
	retV.arr[2] = this->arr[2]/this->arr[3] + v.z;
	retV.arr[3] = 1;
	return retV;
}

Vector Point::operator-(const Point& p)
{
	Vector retV;
	retV.x = this->arr[0]/this->arr[3] - p.arr[0]/p.arr[3];
	retV.y = this->arr[1]/this->arr[3] - p.arr[1]/p.arr[3];
	retV.z = this->arr[2]/this->arr[3] - p.arr[2]/p.arr[3];
	return retV;
}

//code for rotate Vector using Rodrigues' Rotation Formula
void printMatrix(double** matrix);

Vector rotate(Vector x, Vector axis, double angle)
{
    double c = cos((angle/(double)360)*2*pi);
    double s = sin((angle/(double)360)*2*pi);
    Vector retV = (x * c + axis * ((axis,x)*(1-c)) + (axis * x) * s);
    return retV;
}

double** arrToPtr(double arr[4][4])
{
    double** retV = new double*[4];
    for(int i = 0; i<4; i++)
    {
        retV[i] = new double[4];
        for(int j = 0; j<4; j++)
        {
            retV[i][j] = arr[i][j];
        }
    }
    return retV;
}

double** genRotationMat(Vector axis, double angle)
{
    axis.normalize();
    Vector column1 = rotate(*new Vector(1, 0, 0), axis, angle);
    Vector column2 = rotate(*new Vector(0, 1, 0), axis, angle);
    Vector column3 = rotate(*new Vector(0, 0, 1), axis, angle);
    double mat[4][4] ={{column1.x, column2.x, column3.x, 0},
                       {column1.y, column2.y, column3.y, 0},
                       {column1.z, column2.z, column3.z, 0},
                       {0, 0, 0, 1}};
    double** retV = arrToPtr(mat);
    //printMatrix(mat);
    return retV;
}

double** genTranMat(double tx, double ty, double tz)
{
    double tranMat[4][4] = {{1, 0, 0, tx},
                            {0, 1, 0, ty},
                            {0, 0, 1, tz},
                            {0, 0, 0, 1}};
    double** retV = arrToPtr(tranMat);
    return retV;
}

double** genScaleMat(double sx, double sy, double sz)
{
    double scaleMat[4][4] = {{sx, 0, 0, 0},
                             {0, sy, 0, 0},
                             {0, 0, sz, 0},
                             {0, 0, 0, 1}};
    double** retV = arrToPtr(scaleMat);
    return retV;
}

void printMatrix(double** matrix)
{
    cout << "=========================" << endl;
    for(int i = 0; i<4; i++)
    {
        for(int j = 0; j<4; j++)
        {
            cout << fixed << setprecision(3) << matrix[i][j];
            if(j<3) cout << ", ";
        }
        cout << endl;
    }
    cout << "=========================" << endl;
}

double** matMultiply(double** mat1, double** mat2)
{
    double** retV = new double*[4];
    for(int i = 0; i<4; i++)
    {
        retV[i] = new double[4];
        for(int j = 0; j<4; j++)
        {
			retV[i][j] = 0;
            for(int k = 0; k<4; k++)
            {
                retV[i][j]+=mat1[i][k]*mat2[k][j];
            }
        }
    }
    return retV;
}

double** genModelXformMat(Point eye, Point look, Vector up)
{
    Vector l = look - eye;
    l.normalize();
    Vector r = l*up;
    r.normalize();
    Vector u = r*l;
    double rotMat[4][4] = {{r.x, r.y, r.z, 0},
                           {u.x, u.y, u.z, 0},
                           {-l.x, -l.y, -l.z, 0},
                           {0, 0, 0, 1}};
    double** T = genTranMat(-eye.arr[0], -eye.arr[1], -eye.arr[2]);
    double** R = arrToPtr(rotMat);
    //printMatrix(T);
    //printMatrix(R);
    double** V = matMultiply(R,T);
    return V;
}

Point transformPoint(double** TMat, Point p)
{
    Point retV;
    retV.arr[3] = 0;
    for(int i = 0; i<4; i++)
    {
		retV.arr[i] = 0;
        for(int j = 0; j<4; j++)
        {
            retV.arr[i] += TMat[i][j] * p.arr[j];
        }
    }
    for(int i = 0; i<3; i++)
    {
        retV.arr[i] /= retV.arr[3];
    }
    return retV;
}

double** createIDMatrix()
{
    double mat[4][4] = {{1, 0, 0, 0},
                        {0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}};
    double** retV = new double*[4];
    for(int i = 0; i<4; i++)
    {
        retV[i] = new double[4];
        for(int j = 0; j<4; j++)
        {
            retV[i][j] = mat[i][j];
        }
    }
    return retV;
}

class Triangle{
public:
    Point p[3];
    int color[3];
    Triangle operator=(const Triangle& tri)
    {
        for(int i = 0; i<3; i++)
        {
            for(int j = 0; j<4; j++)
            {
                this->p[i].arr[j] = tri.p[i].arr[j];
            }
        }
        return *this;
    }
    bool isYBetween(double ymin, double ymax)
    {
        return !(((p[0].arr[1] < ymin) && (p[1].arr[1] < ymin) && (p[2].arr[1] < ymin)) || ((p[0].arr[1] > ymax) && (p[1].arr[1] > ymax) && (p[2].arr[1] > ymax)));
    }
    bool isXBetween(double xmin, double xmax)
    {
        return !(((p[0].arr[0] < xmin) && (p[1].arr[0] < xmin) && (p[2].arr[0] < xmin)) || ((p[0].arr[0] > xmax) && (p[1].arr[0] > xmax) && (p[2].arr[0] > xmax)));
    }
};

void outputToFile(ofstream& file, Triangle t, int decimalplace)
{
    file << fixed << setprecision(decimalplace) << t.p[0].arr[0] << " " << t.p[0].arr[1] << " " << t.p[0].arr[2] << endl;
    file << fixed << setprecision(decimalplace) << t.p[1].arr[0] << " " << t.p[1].arr[1] << " " << t.p[1].arr[2] << endl;
    file << fixed << setprecision(decimalplace) << t.p[2].arr[0] << " " << t.p[2].arr[1] << " " << t.p[2].arr[2] << endl;
    file << endl;
}
void printTriangle(Triangle t, int decimalplace)
{
//    cout << fixed << setprecision(decimalplace) << t.p[0].arr[0] << " " << t.p[0].arr[1] << " " << t.p[0].arr[2] << endl;
//    cout << fixed << setprecision(decimalplace) << t.p[1].arr[0] << " " << t.p[1].arr[1] << " " << t.p[1].arr[2] << endl;
//    cout << fixed << setprecision(decimalplace) << t.p[2].arr[0] << " " << t.p[2].arr[1] << " " << t.p[2].arr[2] << endl;
//    cout << endl;
      for(int i = 0; i<3; i++)
      {
          t.p[i].print();
      }
}

vector<string> split(string str, char delimiter) {
  vector<string> retV;
  stringstream ss(str);
  string tok;
  while(getline(ss, tok, delimiter)) {
    retV.push_back(tok);
  }
  return retV;
}

bool compare(double pointValue, double distance, int type)
{
    if(type == 1)
    {
        return pointValue <= distance;
    }
    else if(type == 2)
    {
        return pointValue >= distance;
    }
}


vector<Triangle> clipZFar(vector<Triangle> triList, double far)
{
    cout << "here we are" << endl;
    vector<Triangle> retV;
    Point** pIn = new Point*[3];
    for(int i = 0; i<triList.size(); i++)
    {
        //cout << "starting loop" << endl;
        int numOfVertexIn = 0;
        bool b[3];
        b[0] = false;
        b[1] = false;
        b[2] = false;
        if(triList[i].p[0].arr[2] > far )
        {
            //cout << (*triList)[i].p[0].arr[2] << ", " << far << endl;
            numOfVertexIn++;
            b[0] = true;
        }
        if(triList[i].p[1].arr[2] > far )
        {
            //cout << (*triList)[i].p[0].arr[2] << ", " << far << endl;
            numOfVertexIn++;
            b[1] = true;
        }
        if(triList[i].p[2].arr[2] > far )
        {
            //cout << (*triList)[i].p[0].arr[2] << ", " << far << endl;
            numOfVertexIn++;
            b[2] = true;
        }
       // cout << "Points Evaluated" << endl;
        if(numOfVertexIn == 3)
        {
            cout << "what? " << endl;
            Triangle newTri;
            newTri.p[0] = triList[i].p[0];
            newTri.p[1] = triList[i].p[1];
            newTri.p[2] = triList[i].p[2];
            newTri.color[0] = triList[i].color[0];
            newTri.color[1] = triList[i].color[1];
            newTri.color[2] = triList[i].color[2];
            retV.push_back(newTri);
        }
        //cout << "NO PROBLEM!!!!" << endl;
        if(numOfVertexIn == 2)
        {
            if(!b[0])
            {
                pIn[0] = &(triList[i].p[1]);
                pIn[1] = &(triList[i].p[2]);
                pIn[2] = &(triList[i].p[0]);
            }
            if(!b[1])
            {
                pIn[0] = &(triList[i].p[2]);
                pIn[1] = &(triList[i].p[0]);
                pIn[2] = &(triList[i].p[1]);
            }
            if(!b[2])
            {
                pIn[0] = &(triList[i].p[0]);
                pIn[1] = &(triList[i].p[1]);
                pIn[2] = &(triList[i].p[2]);
            }
            double x0 = (*(pIn[0])).getX();
            double x1 = (*(pIn[1])).getX();
            double x2 = (*(pIn[2])).getX();
            double y0 = (*(pIn[0])).getY();
            double y1 = (*(pIn[1])).getY();
            double y2 = (*(pIn[2])).getY();
            double z0 = (*(pIn[0])).getZ();
            double z1 = (*(pIn[1])).getZ();
            double z2 = (*(pIn[2])).getZ();
            Point* p1 = new Point();
            double r1 = (z0 - (far + EPSILON))/(z0 - z2);
            p1->arr[0] = x0 - (x0 - x2) * r1;
            p1->arr[1] = y0 - (y0 - y2) * r1;
            p1->arr[2] = far + EPSILON;
            p1->arr[3] = 1;
           // cout << "p1 found" << endl;
            Point* p2 = new Point();
            double r2 = (z1 - (far + EPSILON))/(z1 - z2);
            p2->arr[0] = x1 - (x1 - x2) * r2;
            p2->arr[1] = y1 - (y1 - y2) * r2;
            p2->arr[2] = far + EPSILON;
            p2->arr[3] = 1;
           // cout << "p2 found" << endl;
            Triangle newTri1;
            newTri1.p[0] = *(pIn[0]);
            newTri1.p[1] = *(pIn[1]);
            newTri1.p[2] = *p1;
            newTri1.color[0] = triList[i].color[0];
            newTri1.color[1] = triList[i].color[1];
            newTri1.color[2] = triList[i].color[2];
            retV.push_back(newTri1);
            Triangle newTri2;
            newTri2.p[0] = *p1;
            newTri2.p[1] = *(pIn[1]);
            newTri2.p[2] = *p2;
            newTri2.color[0] = triList[i].color[0];
            newTri2.color[1] = triList[i].color[1];
            newTri2.color[2] = triList[i].color[2];
            retV.push_back(newTri2);

        }
        //cout << "Still NO PROBLEM!!!!" << endl;
        if(numOfVertexIn == 1)
        {
            cout << "only one is in!" << endl;
            if(b[0])
            {
                pIn[0]=&(triList[i].p[0]);
                pIn[1]=&(triList[i].p[1]);
                pIn[2]=&(triList[i].p[2]);
            }
			if(b[1])
            {
                pIn[0]=&(triList[i].p[1]);
                pIn[1]=&(triList[i].p[0]);
                pIn[2]=&(triList[i].p[2]);
            }
			if(b[2])
            {
                pIn[0]=&(triList[i].p[2]);
                pIn[1]=&(triList[i].p[1]);
                pIn[2]=&(triList[i].p[0]);
            }
            double x0 = (*(pIn[0])).getX();
            double x1 = (*(pIn[1])).getX();
            double x2 = (*(pIn[2])).getX();
            double y0 = (*(pIn[0])).getY();
            double y1 = (*(pIn[1])).getY();
            double y2 = (*(pIn[2])).getY();
            double z0 = (*(pIn[0])).getZ();
            double z1 = (*(pIn[1])).getZ();
            double z2 = (*(pIn[2])).getZ();
            Point* p1 = new Point();
            double r1 = (z0 - (far + EPSILON))/(z0 - z1);
            p1->arr[0] = x0 - (x0 - x1) * r1;
            p1->arr[1] = y0 - (y0 - y1) * r1;
            p1->arr[2] = far + EPSILON;
            p1->arr[3] = 1;
            //cout << "p1 found" << endl;
            Point* p2 = new Point();
            double r2 = (z0 - (far + EPSILON))/(z0 - z2);
            p2->arr[0] = x0 - (x0 - x2) * r2;
            p2->arr[1] = y0 - (y0 - y2) * r2;
            p2->arr[2] = far + EPSILON;
            p2->arr[3] = 1;
            //cout << "p2 found" << endl;
            Triangle newTri;
            newTri.p[0] = *(pIn[0]);
            newTri.p[1] = *p1;
            newTri.p[2] = *p2;
            newTri.color[0] = triList[i].color[0];
            newTri.color[1] = triList[i].color[1];
            newTri.color[2] = triList[i].color[2];
            retV.push_back(newTri);
            //cout << "we have a new triangle here" << endl;

        }

    }
    return retV;

}

vector<Triangle> clipZNear(vector<Triangle> triList, double near)
{
    vector<Triangle> retV;
    Point** pIn = new Point*[3];
    for(int i = 0; i<triList.size(); i++)
    {
        int numOfVertexIn = 0;
        bool b[3] = {false, false, false};
        if(triList[i].p[0].arr[2] < near)
        {
            numOfVertexIn++;
            b[0] = true;
        }
        if(triList[i].p[1].arr[2] < near)
        {
            numOfVertexIn++;
            b[1] = true;
        }
        if(triList[i].p[2].arr[2] < near)
        {
            numOfVertexIn++;
            b[2] = true;
        }
        if(numOfVertexIn == 3)
        {
            Triangle newTri;
            newTri.p[0] = triList[i].p[0];
            newTri.p[1] = triList[i].p[1];
            newTri.p[2] = triList[i].p[2];
            newTri.color[0] = triList[i].color[0];
            newTri.color[1] = triList[i].color[1];
            newTri.color[2] = triList[i].color[2];
            retV.push_back(newTri);
        }
        if(numOfVertexIn == 2)
        {
            if(!b[0])
            {
                pIn[0] = &(triList[i].p[1]);
                pIn[1] = &(triList[i].p[2]);
                pIn[2] = &(triList[i].p[0]);
            }
            if(!b[1])
            {
                pIn[0] = &(triList[i].p[2]);
                pIn[1] = &(triList[i].p[0]);
                pIn[2] = &(triList[i].p[1]);
            }
            if(!b[2])
            {
                pIn[0] = &(triList[i].p[0]);
                pIn[1] = &(triList[i].p[1]);
                pIn[2] = &(triList[i].p[2]);
            }
            double x0 = (*(pIn[0])).getX();
            double x1 = (*(pIn[1])).getX();
            double x2 = (*(pIn[2])).getX();
            double y0 = (*(pIn[0])).getY();
            double y1 = (*(pIn[1])).getY();
            double y2 = (*(pIn[2])).getY();
            double z0 = (*(pIn[0])).getZ();
            double z1 = (*(pIn[1])).getZ();
            double z2 = (*(pIn[2])).getZ();
            Point* p1 = new Point();
            double r1 = (z0 - (near - EPSILON))/(z0 - z2);
            p1->arr[0] = x0 - (x0 - x2) * r1;
            p1->arr[1] = y0 - (y0 - y2) * r1;
            p1->arr[2] = (near - EPSILON);
            p1->arr[3] = 1;
           // cout << "p1 found" << endl;
            Point* p2 = new Point();
            double r2 = (z1 - (near - EPSILON))/(z1 - z2);
            p2->arr[0] = x1 - (x1 - x2) * r2;
            p2->arr[1] = y1 - (y1 - y2) * r2;
            p2->arr[2] = (near - EPSILON);
            p2->arr[3] = 1;
           // cout << "p2 found" << endl;
            Triangle newTri1;
            newTri1.p[0] = *(pIn[0]);
            newTri1.p[1] = *(pIn[1]);
            newTri1.p[2] = *p1;
            newTri1.color[0] = triList[i].color[0];
            newTri1.color[1] = triList[i].color[1];
            newTri1.color[2] = triList[i].color[2];
            retV.push_back(newTri1);
            Triangle newTri2;
            newTri2.p[0] = *p1;
            newTri2.p[1] = *(pIn[1]);
            newTri2.p[2] = *p2;
            newTri2.color[0] = triList[i].color[0];
            newTri2.color[1] = triList[i].color[1];
            newTri2.color[2] = triList[i].color[2];
            retV.push_back(newTri2);

        }
        if(numOfVertexIn == 1)
        {
            if(b[0])
            {
                pIn[0]=&(triList[i].p[0]);
                pIn[1]=&(triList[i].p[1]);
                pIn[2]=&(triList[i].p[2]);
            }
			if(b[1])
            {
                pIn[0]=&(triList[i].p[1]);
                pIn[1]=&(triList[i].p[0]);
                pIn[2]=&(triList[i].p[2]);
            }
			if(b[2])
            {
                pIn[0]=&(triList[i].p[2]);
                pIn[1]=&(triList[i].p[1]);
                pIn[2]=&(triList[i].p[0]);
            }
            double x0 = (*(pIn[0])).getX();
            double x1 = (*(pIn[1])).getX();
            double x2 = (*(pIn[2])).getX();
            double y0 = (*(pIn[0])).getY();
            double y1 = (*(pIn[1])).getY();
            double y2 = (*(pIn[2])).getY();
            double z0 = (*(pIn[0])).getZ();
            double z1 = (*(pIn[1])).getZ();
            double z2 = (*(pIn[2])).getZ();
            Point* p1 = new Point();
            double r1 = (z0 - (near - EPSILON))/(z0 - z1);
            p1->arr[0] = x0 - (x0 - x1) * r1;
            p1->arr[1] = y0 - (y0 - y1) * r1;
            p1->arr[2] = (near - EPSILON);
            p1->arr[3] = 1;
           // cout << "p1 found" << endl;
            Point* p2 = new Point();
            double r2 = (z0 - (near - EPSILON))/(z0 - z2);
            p2->arr[0] = x0 - (x0 - x2) * r2;
            p2->arr[1] = y0 - (y0 - y2) * r2;
            p2->arr[2] = (near - EPSILON);
            p2->arr[3] = 1;
           // cout << "p2 found" << endl;
            Triangle newTri;
            newTri.p[0] = *(pIn[0]);
            newTri.p[1] = *p1;
            newTri.p[2] = *p2;
            newTri.color[0] = triList[i].color[0];
            newTri.color[1] = triList[i].color[1];
            newTri.color[2] = triList[i].color[2];
            retV.push_back(newTri);

        }

    }
    return retV;

}


vector<Triangle> clipZ(vector<Triangle> triList, double near, double far)
{
    return clipZFar(clipZNear(triList, near), far);
}

vector<Point> findIntersects(Triangle tri, double xline)
{
    vector<Point> retV;
    double t1 = (tri.p[0].arr[0] - xline) / (tri.p[0].arr[0] - tri.p[1].arr[0]); //01
    double t2 = (tri.p[1].arr[0] - xline) / (tri.p[1].arr[0] - tri.p[2].arr[0]); //12
    double t3 = (tri.p[2].arr[0] - xline) / (tri.p[2].arr[0] - tri.p[0].arr[0]); //20
    Point P1;
    P1.arr[0] = xline;
    P1.arr[1] = tri.p[0].arr[1] - (tri.p[0].arr[1] - tri.p[1].arr[1]) * t1;
    P1.arr[2] = tri.p[0].arr[2] - (tri.p[0].arr[2] - tri.p[1].arr[2]) * t1;
    P1.arr[3] = 1;
    Point P2;
    P2.arr[0] = xline;
    P2.arr[1] = tri.p[1].arr[1] - (tri.p[1].arr[1] - tri.p[2].arr[1]) * t2;
    P2.arr[2] = tri.p[1].arr[2] - (tri.p[1].arr[2] - tri.p[2].arr[2]) * t2;
    P2.arr[3] = 1;
    Point P3;
    P3.arr[0] = xline;
    P3.arr[1] = tri.p[2].arr[1] - (tri.p[2].arr[1] - tri.p[0].arr[1]) * t3;
    P3.arr[2] = tri.p[2].arr[2] - (tri.p[2].arr[2] - tri.p[0].arr[2]) * t3;
    P3.arr[3] = 1;
    if(0 <= t1 && t1 <= 1)
    {
        retV.push_back(P1);
    }
    if(0 <= t2 && t2 <= 1)
    {
        retV.push_back(P2);
    }
    if(0 <= t3 && t3 <= 1)
    {
        retV.push_back(P3);
    }
    if(retV.size() == 1)
    {
        cout << "(t1, t2, t3) = (" << t1 << ", " << t2 << ", " << t3 << ")" << endl;
    }
    if(retV.size() == 3)
    {
        cout << "What's happening?" << endl;
    }
    return retV;
}

void printTable(double** table, int tableRows, int tableCols)
{
    for(int i = 0; i<tableRows; i++)
    {
        for(int j = 0; j<tableCols; j++)
        {
            cout << trunc(table[i][j]) << " ";
        }
        cout << endl;
    }
}

double** makeZBuffer(Triangle t, int screenWidth, int screenHeight)
{
    double** retV = new double*[screenWidth];
    for(int i = 0; i<screenWidth; i++)
    {
        retV[i] = new double[screenHeight];
        for(int j = 0; j<screenHeight; j++)
        {
            retV[i][j] = 1.0;
        }
    }
    for(int i = 0; i<screenWidth; i++)
    {
        double xline = -1.0 + (double)(2*i + 1)/(double)screenWidth;
       // cout << "finding intersect for i = " << i << endl;
        vector<Point> intsect = findIntersects(t, xline);
       // cout << "intersect found" << endl;
        for(int j = 0; j<screenHeight; j++)
        {
            if(intsect.size() == 2)
            {

                double yline = -1.0 +(double)(2*j + 1)/(double)screenHeight;
                if(intsect[1].arr[1] > intsect[0].arr[1])
                {
                    double c = ((intsect[1].arr[1] - yline)/(intsect[1].arr[1] - intsect[0].arr[1]));
                    if(0 <= c && c <= 1) retV[i][screenHeight - j] = intsect[1].arr[2] - (intsect[1].arr[2] - intsect[0].arr[2])*c;
                }
                else
                {
                    double c = ((intsect[0].arr[1] - yline)/(intsect[0].arr[1] - intsect[1].arr[1]));
                    if(0 <= c && c <= 1) retV[i][screenHeight - j] = intsect[0].arr[2] - (intsect[0].arr[2] - intsect[1].arr[2])*c;
                }
            }
        }
    }
    return retV;

}
void writeImage(vector<Triangle>* triList, int screenWidth, int screenHeight, int* backColor)
{
    bitmap_image image(screenWidth, screenHeight);
    double*** zBuffers = new double**[triList->size()];
   // cout << "here" << endl;
    double** globalBuffer = new double*[screenWidth];
    for(int i = 0; i<screenWidth; i++)
    {
        globalBuffer[i] = new double[screenHeight];
        for(int j = 0; j<screenHeight; j++)
        {
            globalBuffer[i][j] = 1.0;
            image.set_pixel(i, j, backColor[0], backColor[1], backColor[2]);
        }
    }
   // cout << "oh yeah!" << endl;
    int lsize = triList->size();
    cout << "What is the size? " << lsize << endl;
    for(int t = 0; t < lsize ; t++)
    {
        //cout << " making t  = " << t << endl;
        if(!(*triList)[t].isYBetween(-1, 1) || !(*triList)[t].isXBetween(-1, 1))
        {
            continue;
        }
        zBuffers[t] = makeZBuffer((*triList)[t], screenWidth, screenHeight);
        //cout << " t  = " << t << " made" << endl;
        for(int i = 0; i<screenWidth; i++)
        {
            for(int j = 0; j<screenHeight; j++)
            {
                if(zBuffers[t][i][j] > -1.0)
                {
                    if(zBuffers[t][i][j] < globalBuffer[i][j])
                    {
                        image.set_pixel(i, j, (*triList)[t].color[0], (*triList)[t].color[1], (*triList)[t].color[2]);
                        globalBuffer[i][j] = zBuffers[t][i][j];
                    }
                }
            }
        }
    }
    image.save_image("out.bmp");

}
/* void destroyMat(double** mat)
{
    for(int i = 0; i<4; i++)
    {
        delete[] mat[i];
    }
    delete[] mat;
} */

void printTriList(vector<Triangle> triList)
{
    for(int i = 0; i<triList.size(); i++)
    {
        cout << endl;
        printTriangle(triList[i], 7);
        cout << endl;
    }
}

int main()
{
   //double** test = genRotationMat(*new Vector(1,0,1), 90);
   //printMatrix(test);
//    double** test2 = genModelXformMat(*new Point(0, 0, 50, 1), *new Point(), *new Vector(0, 1, 0));
//    printMatrix(test2);
    Point look;
    Point eye;
    Vector up;
    ifstream scene;
    ofstream stage1;
    ofstream stage2;
    ofstream stage3;
    scene.open("scene.txt");
    stage1.open("stage1.txt");
    string line;
    vector<string> elems;
    vector<string> gluPerspective;
    int linecount = 0;
    stack<double**> stk;
    vector<Triangle> triList;
    Triangle* tri;
    int triangle = 4;
    bool transProcess = false;
    bool rotProcess = false;
    bool scaleProcess = false;
    int screenHeight = 0;
    int screenWidth = 0;
    int backgroundColor[3];
//    vector<double**> forDelete;
    stk.push(createIDMatrix());
    while(true)
    {
        getline(scene, line);
        linecount++;
        elems = split(line, ' ');
        if(linecount == 1)
        {
            eye.arr[0] = atof(elems[0].c_str());
            eye.arr[1] = atof(elems[1].c_str());
            eye.arr[2] = atof(elems[2].c_str());
            eye.arr[3] = 1;
        }
        else if(linecount == 2)
        {
            look.arr[0] = atof(elems[0].c_str());
            look.arr[1] = atof(elems[1].c_str());
            look.arr[2] = atof(elems[2].c_str());
            look.arr[3] = 1;
        }
        else if(linecount == 3)
        {
            up.x = atof(elems[0].c_str());
            up.y = atof(elems[1].c_str());
            up.z = atof(elems[2].c_str());
        }
        else if(linecount == 4)
        {
            gluPerspective = elems;
        }
        else if(linecount == 5)
        {
            screenWidth = atoi(elems[0].c_str());
            screenHeight = atoi(elems[1].c_str());
            //cout << screenWidth << " " << screenHeight << endl;
        }
        else if(linecount == 6)
        {
            for(int i = 0; i<3; i++)
            {
                backgroundColor[i] = atoi(elems[i].c_str());
                //cout << backgroundColor[i] << " ";
            }
            //cout << endl;
        }
        else
        {
            if(triangle < 3)
            {
                Point p(atof(elems[0].c_str()), atof(elems[1].c_str()), atof(elems[2].c_str()), 1);
                tri->p[triangle] = transformPoint(stk.top(), p);
                triangle++;
            }
            else if(triangle == 3)
            {
                for(int i = 0; i<3; i++)
                {
                    tri->color[i] = atoi(elems[i].c_str());
                }
                outputToFile(stage1, *tri, 7);
                triList.push_back(*tri);
                triangle++;
            }
            else if(line.compare("triangle") == 0)
            {
                tri = new Triangle();
                triangle = 0;
            }
            else if(line.compare("translate") == 0)
            {
                transProcess = true;
            }
            else if(line.compare("scale") == 0)
            {
                scaleProcess = true;
            }
            else if(line.compare("rotate") == 0)
            {
                rotProcess = true;
            }
            else if(line.compare("push") == 0)
            {
                double** top = stk.top();
                double topC[4][4];
                for(int i = 0; i<4; i++)
                {
                    for(int j = 0; j<4; j++)
                    {
                        topC[i][j] = top[i][j];
                    }
                }
                stk.push(arrToPtr(topC));
            }
            else if(line.compare("pop") == 0)
            {
                stk.pop();
            }
            else if(line.compare("end") == 0)
            {
                break;
            }
            else if(transProcess == true)
            {
                double** tMat = genTranMat(atof(elems[0].c_str()), atof(elems[1].c_str()), atof(elems[2].c_str()));
                //cout << line << endl;
                //printMatrix(tMat);
                double** stkTopCopy = stk.top();
                stk.top() = matMultiply(stkTopCopy, tMat);
//                forDelete.push_back(stkTopCopy);
//                forDelete.push_back(tMat);
                transProcess = false;
            }
            else if(scaleProcess == true)
            {
                double** sMat = genScaleMat(atof(elems[0].c_str()), atof(elems[1].c_str()), atof(elems[2].c_str()));
                //cout << line << endl;
                //printMatrix(sMat);
                double** stkTopCopy = stk.top();
                stk.top() = matMultiply(stkTopCopy, sMat);
//                forDelete.push_back(stkTopCopy);
//                forDelete.push_back(sMat);
                scaleProcess = false;
            }
            else if(rotProcess == true)
            {
                double angle = atof(elems[0].c_str());
                Vector axis(atof(elems[1].c_str()), atof(elems[2].c_str()), atof(elems[3].c_str()));
                double** rMat = genRotationMat(axis, angle);
                //cout << line << endl;
                //printMatrix(rMat);
                double** stkTopCopy = stk.top();
                stk.top() = matMultiply(stkTopCopy, rMat);
//                forDelete.push_back(stkTopCopy);
//                forDelete.push_back(rMat);
                rotProcess = false;
            }
        }
    }
    stage1.close();
    scene.close();
    stage2.open("stage2.txt");
    double** modelXForm = genModelXformMat(eye,look,up);
    //forDelete.push_back(modelXForm);
    for(int i = 0; i<triList.size(); i++)
    {
        for(int j = 0; j<3; j++)
        {
            triList[i].p[j] = transformPoint(modelXForm, triList[i].p[j]);
        }
        outputToFile(stage2, triList[i], 7);
    }
    stage2.close();
//    forDelete.push_back(modelXForm);
    stage3.open("stage3.txt");
    double fovY = atof(gluPerspective[0].c_str());
    double aspectRatio = atof(gluPerspective[1].c_str());
    double near = atof(gluPerspective[2].c_str());
    double far = atof(gluPerspective[3].c_str());
    triList = clipZ(triList, -near, -far);
    cout << "LOOK HERE, triList size: " << triList.size() << endl;
    double fovX = fovY * aspectRatio;
    cout << fovY << " " << aspectRatio << endl;
    cout << fovX << endl;
    double t = near * tan(((fovY/(double)(2*360)))*2*pi);
    double r = near * tan(((fovX/(double)(2*360)))*2*pi);
    cout << t << " " << r << endl;
    double P_arr[4][4] = {{near/r, 0, 0, 0},
                          {0, near/t, 0, 0},
                          {0, 0, -(far+near)/(far-near), -(2*far*near)/(far-near)},
                          {0, 0, -1, 0}};
    double** P = arrToPtr(P_arr);
    for(int i = 0; i<triList.size(); i++)
    {
        for(int j = 0; j<3; j++)
        {
            triList[i].p[j] = transformPoint(P, triList[i].p[j]);
        }
        outputToFile(stage3, triList[i], 7);
    }
    writeImage(&triList, screenWidth, screenHeight, backgroundColor);
    stage3.close();
    return 0;
}
