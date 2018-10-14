//
//  Class_Tensor.hpp
//  Proj_EffectiveDiffusionNEW
//
//  Created by Raffaele Marino on 19/10/15.
//  Copyright (c) 2015 Raffaele Marino. All rights reserved.
//

#ifndef Proj_EffectiveDiffusion_Class_Tensor_hpp
#define Proj_EffectiveDiffusion_Class_Tensor_hpp
#include "header_s.h"

class Tensor{
    
public:
	
     Tensor(){};
    ~Tensor(){};
    vector<double> NewTensorBodysTRANS( char * argv[], int &dim);
    vector<double> NewTensorBodysROT(char * argv[], int &dim);
    void Display();
    double ElementTensor(int i, int j);
    double Trace();
    double  operator [] (int i);
    double InvTrace();
    vector<double> InvVec();
    
    
    
private:
    
    int _dim;
    int _myi;
    double _tracegammameno1;
    vector<vector<double> > _vec;
    vector<double> _myvec;
    vector<double> _myvecsqrt;
   void _SqrtMatrix();
    
};


inline double Tensor::Trace(){return (_vec[0][0]+_vec[1][1]+_vec[2][2]);};

inline double Tensor::InvTrace(){return _tracegammameno1;};

inline vector<double> Tensor::InvVec(){return _myvecsqrt;};
 
inline double Tensor::ElementTensor(int i, int j){ return _vec[i][j];};

vector<double> Tensor::NewTensorBodysTRANS(char * argv[], int &dim){
    vector<double> vecelemdiag;
    vector<double> vecelemoutdiag;
    _dim=dim;
    _vec.resize(_dim);
    for (int i=0; i<_dim; ++i) {
        _vec[i].resize(_dim);
    }
    _myi=2;    /*Value to use when we have the initial value from external INPUT */
    

    vecelemdiag.push_back(stod(argv[3]));
    vecelemdiag.push_back(stod(argv[6]));
    vecelemdiag.push_back(stod(argv[7]));
    vecelemdiag.push_back(stod(argv[6]));
    vecelemdiag.push_back(stod(argv[4]));
    vecelemdiag.push_back(stod(argv[8]));
    vecelemdiag.push_back(stod(argv[7]));
    vecelemdiag.push_back(stod(argv[8]));
    vecelemdiag.push_back(stod(argv[5]));
    
    for (int i=0; i<_vec.size(); ++i) { /*   This function can be optimized, but we can do that later*/
        for (int j=0; j<_vec[i].size(); ++j) {
          
                _vec[i][j]=vecelemdiag[3*i+j]; //This is given in micrometer
                
            cout<<_vec[i][j]<<" ";
        
        }
        cout<<endl;
    }
  
    _SqrtMatrix();
    
    return _myvec;
};


vector<double> Tensor::NewTensorBodysROT(char * argv[], int &dim){
    vector<double> vecelemdiag;
    vector<double> vecelemoutdiag;
    _dim=dim;
    _vec.resize(_dim);
    for (int i=0; i<_dim; ++i) {
        _vec[i].resize(_dim);
    }
    _myi=8;
 
    for (int i=1; i<_dim+1; ++i) {
        vecelemdiag.push_back(stod(argv[_myi+i]));
    }
    
    for (int i=0; i<_vec.size(); ++i) {  /*   This function can be optimized, but we can do that later*/
        for (int j=0; j<_vec[i].size(); ++j) {
            if (i==j && vecelemdiag[i]!=0. ) {
	          _vec[i][j]=1./vecelemdiag[i];
	          cout<<_vec[i][j]<<endl;
            }
            _myvec.push_back(_vec[i][j]);
        }
    }
    return _myvec;
};

/* In the last two function written above we do not loose a lot of time in the run time.*/

inline double  Tensor:: operator [] (int i){
    /*if(i==_myvec.size()){ //fast version
      cout<<"out of size";
        return _myvec[0];
    }*/
    return  _myvec[i];
};

void Tensor::Display(){
    for (int i=0; i<_myvec.size(); ++i) {
        cout<<_myvec[i]<<endl;
    }
};



void Tensor::_SqrtMatrix(){
    mat A(_dim,_dim);
    mat R(_dim,_dim);
    
            A << _vec[0][0]<< _vec[0][1]<< _vec[0][2] << endr
            << _vec[1][0]<< _vec[1][1]<< _vec[1][2] << endr
            << _vec[2][0]<< _vec[2][1]<< _vec[2][2] << endr;
            

    
    A.print("A:");
    R=inv(A);
    cout << "det(A): " << det(A) << endl;
    cout << "inv(A): " << endl << R << endl;
    A.save("A.txt", raw_ascii);
    R.save("R.txt", raw_ascii);
    vec eigval;
    mat eigvec;
    eig_sym( eigval, eigvec, R );
    eigval.print("eigval:");
    eigvec.print("eignvec:");
    mat B=eigvec;
    B.print("B:");
    mat I=inv(B);
    I.print("I:");
    mat D=I*R;
    D=D*B;
    D.print("D:");
    mat C=sqrt(D);
    C.print("C:");
    mat S=B*C;
    S=S*I;
    S.print("S:");
    S.save("S.txt", raw_ascii);
    mat N=S*S;
    N.print("N:");
    mat O=R-N;
    for (int i=0; i<_dim; ++i) {
        for (int j=0; j<_dim; ++j) {
            if(abs(O[i*_dim+j])>0.0000001){
                cout<<"BAD MATRIX! Change the values of Translational Friction Tensor. If the issue continues you can contact the author: raffaele_marino@nordita.org"<<endl;
                exit(-1);}
            _myvecsqrt.push_back(S[i*_dim+j]);
            _myvec.push_back(R[i*_dim+j]);
        }
    }
    cout<<"Nice Matrix! The program starts to integer"<<endl;
    _tracegammameno1=trace(R);
    
};

#endif
