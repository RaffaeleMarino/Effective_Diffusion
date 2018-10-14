//
//  Class_Quaternion.hpp
//  Proj_EffectiveDiffusion
//
//  Created by Raffaele Marino on 19/10/15.
//  Copyright (c) 2015 Raffaele Marino. All rights reserved.
//


#ifndef Proj_EffectiveDiffusion_Class_Quaternion_hpp
#define Proj_EffectiveDiffusion_Class_Quaternion_hpp
#include "header_s.h"

class Quaternion{
public:
	
    Quaternion(void)
    {
        _x = 0;
        _y = 0;
        _z = 0;
        _w = 0;
    }
    
    Quaternion(double wi, double xi, double yi, double zi)
    {
        _w = wi;
        _x = xi;
        _y = yi;
        _z = zi;
    }

    Quaternion(double v[4])
    {
        _w = v[0];
        _x = v[1];
        _y = v[2];
        _z = v[3];
    }

    Quaternion(const Quaternion& q)
    {
        _w = q._w;
        _x = q._x;
        _y = q._y;
        _z = q._z;
    }
    
    ~Quaternion(){}
    
    void Prodforonehalf(double s);
    void RotationQuaternion();
    void NormalizationQuaternion();
    void PushQuaternion(double &a, vector<double> &b);
    void Display();
    bool operator != (const Quaternion& q);
    bool operator == (const Quaternion& q);
    double norm();
    double magnitude();
    double ScalarQuaternion();
    double xQuaternion();
    double yQuaternion();
    double zQuaternion();
    vector<double> VectorQuaternion();
    Quaternion scale( double s);
    Quaternion inverse();
    Quaternion conjugate();
    Quaternion operator = (const Quaternion& q);
    Quaternion operator + (const Quaternion& q);
    Quaternion operator - (const Quaternion& q);
    Quaternion operator * (const Quaternion& q);
    Quaternion operator / (Quaternion& q);
    Quaternion&  operator += (const Quaternion& q);
    Quaternion&  operator -= (const Quaternion& q);
    Quaternion&  operator *= (const Quaternion& q);
    Quaternion&  operator /= (Quaternion& q);
    
private:
    
   	double _w, _x, _y, _z;
    
    vector<double> _b;
	
	
};


inline Quaternion Quaternion::operator = (const Quaternion& q)
{
    _w = q._w;
    _x = q._x;
    _y = q._y;
    _z = q._z;
    
    return (*this);
};

inline Quaternion Quaternion::operator + (const Quaternion& q)
{
    return Quaternion(_w+q._w, _x+q._x, _y+q._y, _z+q._z);
};

inline Quaternion Quaternion::operator - (const Quaternion& q)
{
    return Quaternion(_w-q._w, _x-q._x, _y-q._y, _z-q._z);
};

 inline Quaternion Quaternion::operator * (const Quaternion& q)
{
    return Quaternion(
                      _w*q._w - _x*q._x - _y*q._y - _z*q._z,
                      _w*q._x + _x*q._w + _y*q._z - _z*q._y,
                      _w*q._y + _y*q._w + _z*q._x - _x*q._z,
                      _w*q._z + _z*q._w + _x*q._y - _y*q._x);
};

 inline Quaternion Quaternion::operator / (Quaternion& q)
{
    return ((*this) * (q.inverse()));
};


Quaternion& Quaternion::operator += (const Quaternion& q)
{
    _w += q._w;
    _x += q._x;
    _y += q._y;
    _z += q._z;
    
    return (*this);
};


Quaternion& Quaternion::operator -= (const Quaternion& q)
{
    _w -= q._w;
    _x -= q._x;
    _y -= q._y;
    _z -= q._z;
    
    return (*this);
};


Quaternion& Quaternion::operator *= (const Quaternion& q)
{
    double w_val = _w*q._w - _x*q._x - _y*q._y - _z*q._z;
    double x_val = _w*q._x + _x*q._w + _y*q._z - _z*q._y;
    double y_val = _w*q._y + _y*q._w + _z*q._x - _x*q._z;
    double z_val = _w*q._z + _z*q._w + _x*q._y - _y*q._x;
    
    _w = w_val;
    _x = x_val;
    _y = y_val;
    _z = z_val;
    
    return (*this);
};


Quaternion& Quaternion::operator /= (Quaternion& q)
{
    (*this) = (*this)*q.inverse();
    return (*this);
};

inline bool Quaternion::operator != (const Quaternion& q)
{
    return (_w!=q._w || _x!=q._x || _y!=q._y || _z!=q._z) ? true : false;
};


inline bool Quaternion::operator == (const Quaternion& q)
{
    return (_w==q._w && _x==q._x && _y==q._y && _z==q._z) ? true : false;
};


double Quaternion::norm()
{
    return (_w*_w + _x*_x + _y*_y + _z*_z);
};


double Quaternion::magnitude()
{
    return sqrt(norm());
};

 //Using for moltipling the scalar
Quaternion Quaternion::scale(double  s)
{
    return Quaternion(_w*s, _x*s, _y*s, _z*s);
};

 //Using for moltipling the scalar
void Quaternion::Prodforonehalf(double  s)
{
    _w=_w*s;
    _x=_x*s;
    _y=_y*s;
    _z=_z*s;
};


Quaternion Quaternion::inverse()
{
    return conjugate().scale(1./norm());
};


Quaternion Quaternion::conjugate()
{
    return Quaternion(_w, -_x, -_y, -_z);
};


void Quaternion::NormalizationQuaternion(){
    double _magnitude=magnitude();
    _w=_w/_magnitude;
    _x=_x/_magnitude;
    _y=_y/_magnitude;
    _z=_z/_magnitude;
    
};


void Quaternion::PushQuaternion(double &a, vector<double> &b){
    _w=a;
    _x=b[0];
    _y=b[1];
    _z=b[2];
};

void Quaternion::RotationQuaternion(){

    double teta;
    teta=PI/4.; // See what Ralf needs for the paper.
    _w=cos(teta/2.);
    _x=0.;
    _y=0.;
    _z=sin(teta/2.);
    
};

 inline double Quaternion::ScalarQuaternion(){return _w;};

 inline double Quaternion::xQuaternion(){return _x;};

 inline double Quaternion::yQuaternion(){return _y;};

 inline double Quaternion::zQuaternion(){return _z;};



vector<double> Quaternion::VectorQuaternion(){
    _b.resize(3);
    _b[0]=_x;
    _b[1]=_y;
    _b[2]=_z;
    return _b;};


void Quaternion::Display(){
    cout<<_w<<" "<<_x<<" "<<_y<<" "<<_z<<endl;
};

#endif
