//
//  Class_Temperature.h
//  Proj_EffectiveDiffusion
//
//  Created by Raffaele Marino on 15/06/15.
//  Copyright (c) 2015 Raffaele Marino. All rights reserved.
//

#ifndef Proj_EffectiveDiffusion_Class_Temperature_hpp
#define Proj_EffectiveDiffusion_Class_Temperature_hpp
#include "header_s.h"

class Temperature{
    
public:
	
    Temperature(){};
    ~Temperature(){};
    const double ConstTemperature(); /*This function returns the value of constant temperature that we give as INPUT */
    void Inizialization(char * argv[],int &dim);
    vector<double> GradientTemperature(vector<double> &a);
    
    
    
    
private:
    int _dim;
    vector<double> _vec;
    double _Temp;
    double _FormulaGrad(double &a);
};

inline void Temperature::Inizialization(char * argv[],int &dim){
    _dim=dim;
    _vec.resize(_dim);
    int k=0;
    k=18;
    _Temp=stod(argv[k]);
}

inline const double Temperature::ConstTemperature(){
    double Temp;
    Temp=_Temp;  /*Kelvin's scale*/
    return Temp;
};

/*These following function has to debugged if we want to use a force dipendent by the position*/
double Temperature::_FormulaGrad(double &a){
    return ConstTemperature()*a;
};

inline vector<double> Temperature::GradientTemperature(vector<double> &a){
    double _a;
    for (int i=0; i<_dim; ++i) {
        _a=a[i];
        _vec[i]=_FormulaGrad(_a);
    }
    
    return _vec;
    
}

#endif
