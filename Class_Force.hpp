//
//  Class_Force.hpp
//  Proj_EffectiveDiffusionNEW
//
//  Created by Raffaele Marino on 19/10/15.
//  Copyright (c) 2015 Raffaele Marino. All rights reserved.
//

#ifndef Proj_EffectiveDiffusion_Class_Force_hpp
#define Proj_EffectiveDiffusion_Class_Force_hpp
#include "header_s.h"

class Force{
    
public:
	
    Force(){};
    ~Force(){};
    const vector<double> ConstForce();
    void Inizialization(char * argv[],int &dim);
    void Display();
    vector<double> GradientForce(vector<double> &a); //This function has to be implemented.
    double  operator [] (int i);
    
    
    
    
private:
    
    int _dim;
    vector<double> _vec;
    double _FormulaGrad( int i, double &a);

};

inline void Force::Inizialization(char * argv[],int &dim){
    _dim=dim;
    _vec.resize(_dim);
    int k=0;
    k=12;
    for (int i=0; i<_dim; ++i) {
      _vec[i]=stod(argv[k+i]);
    }
};

inline double  Force:: operator [] (int i){
    if(i==_vec.size()){
        cout<<"out of size";
        exit(-1);
        return _vec[0];
    }
    return  _vec[i];
};

inline const vector<double> Force::ConstForce(){ return _vec;};


double Force::_FormulaGrad(int i,double &a){
    return _vec[i]*a;
};

vector<double> Force::GradientForce(vector<double> &a){
    double _a;
    for (int i=0; i<_dim; ++i) {
        _a=a[i];
        _vec[i]=_FormulaGrad(i,_a);
    }
    
    return _vec;
    
}

void Force::Display(){
    for (int i=0; i<_dim; ++i) {
        cout<<_vec[i]<<endl;
    }
};
#endif
