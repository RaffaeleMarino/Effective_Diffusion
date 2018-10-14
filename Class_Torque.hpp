//
//  Class_Torque.hpp
//  Proj_EffectiveDiffusion
//
//  Created by Raffaele Marino on 19/10/15.
//  Copyright (c) 2015 Raffaele Marino. All rights reserved.
//

#ifndef Proj_EffectiveDiffusion_Class_Torque_hpp
#define Proj_EffectiveDiffusion_Class_Torque_hpp
#include "header_s.h"

class Torque{
    
public:
    
    Torque(){};
    ~Torque(){};
    const vector<double> ConstTorque();
    void Inizialization(char * argv[],int &dim);
    void Display();
    vector<double> GradientTorque(vector<double> &a); //This functionn has to be implemented.
    double  operator [] (int i);
    
    
    
private:

    int _dim;
    vector<double> _vec;
    double _FormulaGrad( int i, double &a);
    
};

inline void Torque::Inizialization(char * argv[],int &dim){
    _dim=dim;
    _vec.resize(_dim);
    int k=0;
    k=15;
    for (int i=0; i<_dim; ++i) {
        _vec[i]=stod(argv[k+i]);
    }
};

inline double  Torque:: operator [] (int i){
    if(i==_vec.size()){
        cout<<"out of size";
	exit(-1);
        return _vec[0];
    }
    return  _vec[i];
};

inline const vector<double> Torque::ConstTorque(){ return _vec;};

/*These following two function has to be improved and debugged when we want to use torque dipendent by the position*/
double Torque::_FormulaGrad(int i,double &a){
    return _vec[i]*a;
};

inline vector<double> Torque::GradientTorque(vector<double> &a){
    double _a;
    for (int i=0; i<_dim; ++i) {
        _a=a[i];
        _vec[i]=_FormulaGrad(i, _a);
    }
    return _vec;  
};

void Torque::Display(){
    for (int i=0; i<_dim; ++i) {
        cout<<_vec[i]<<endl;
    }
};
#endif
