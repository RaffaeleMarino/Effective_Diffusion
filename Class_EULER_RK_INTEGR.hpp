//
//  Class_EULER_RK_INTEGR.h
//  Proj_EffectiveDiffusionNEW
//
//  Created by Raffaele Marino on 19/10/15.
//  Copyright (c) 2015 Raffaele Marino. All rights reserved.
//

#ifndef Proj_EffectiveDiffusion_Class_EULER_RK_INTEGR_hpp
#define Proj_EffectiveDiffusion_Class_EULER_RK_INTEGR_hpp
#include "header_s.h"
#include "Class_Temperature.hpp"
#include "Class_Tensor.hpp"
#include "Class_Quaternion.hpp"
#include "Class_Force.hpp"
#include "Class_Torque.hpp"





class Eulero{
	
public:
	
	Eulero(int &dim){
        _dim=dim;
        _vectorposition.resize(_dim);
        _col=_dim;
        
    };
    ~Eulero(){};
  void Integration (vector<double> & posizionextempot, vector<double> & posizioneytempot, vector<double> & posizioneztempot, Tensor &Thetensorgamma, Tensor &Thetensoreta, Temperature &TheTemperature, Force &TheForce, Torque &TheTorque);
    void PrintonFile ();
    void Mean (vector<double> & xvec,vector<double> &yvec,vector<double> &zvec);
    
private:
    
    vector<double> _vectorposition;
    vector<double> _vecmeanpassix;
    vector<double> _vecmeanpassiy;
    vector<double> _vecmeanpassiz;
    vector<double> _vecmeanpassixsquare;
    vector<double> _vecmeanpassiysquare;
    vector<double> _vecmeanpassizsquare;
    vector<double> _vecmeanpassixysquare;
    vector<double> _vecmeanpassiyzsquare;
    vector<double> _vecmeanpassizxsquare;
    vector<double> _vecmeanpassixxx;
    vector<double> _vecmeanpassiyyy;
    vector<double> _vecmeanpassizzz;

      vector<double> _vecmeanpassixxxx;
    vector<double> _vecmeanpassiyyyy;
    vector<double> _vecmeanpassizzzz;
  
    vector<double> _vecmeanpassixyx;
    vector<double> _vecmeanpassixyz;
    vector<double> _vecmeanpassixzx;
     vector<double> _vecmeanpassizzy;
    vector<double> _vecmeanpassiyxy;
    vector<double> _vecmeanpassiyxz;
    vector<double> _vecmeanpassiyzy;
    vector<double> _vecmeanpassizxz;
    vector<double> _vecmeanpassizxy;
    double _constT;
    int _col;
    int _dim;
  string _mystring;
  ostringstream _mystringgammax, _mystringgammay, _mystringgammaz, _mystringforcex, _mystringforcey, _mystringforcez, _mystringgammaxy, _mystringgammaxz,_mystringgammayz;
    
    unsigned int _GetRandom();
   void PrintonFileTrajectory(vector<double> &xvec,vector<double> &yvec,vector<double> &zvec);
  vector<double> KeffRS(Tensor &Thetensorgamma, Tensor &Thetensoreta, Temperature &TheTemperature, Force &TheForce);
    vector<double> _Random_Numb_Gen();
    
    
};

vector<double> Eulero::KeffRS(Tensor &Thetensorgamma, Tensor &Thetensoreta, Temperature &TheTemperature, Force &TheForce){
  //here everything is dyagonal
  
  const double constT=TheTemperature.ConstTemperature();
  const double traceeta=(1./Thetensoreta[0])+(1./Thetensoreta[4])+(1./Thetensoreta[8]);
  const double constduekT=1./(double)(2.*k_B*constT);
const double constduebkT=1./(double)(2.*k_B*constT*(Thetensoreta[8]*Thetensoreta[0]+Thetensoreta[0]*Thetensoreta[4]+Thetensoreta[4]*Thetensoreta[8]));
  const double tracegammam1=Thetensorgamma.InvTrace();
const double invtraceeta=Thetensoreta.Trace();
  vector<double> vecforce;
  vector<double> veckeff;
  vecforce=TheForce.ConstForce();
    ostringstream T_str;
    T_str<<_constT;
    string mystring;
    const string title4="KeffRS";
    const string dat4=".txt";
    mystring=title4;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();
    mystring+=T_str.str();
    mystring+=dat4;
    ofstream outfilew7(mystring.c_str(), ios_base::app );
  

    const double g=constduebkT*((Thetensorgamma[0]-(tracegammam1/3.))*(Thetensorgamma[0]-(tracegammam1/3.))*Thetensoreta[0]+(Thetensorgamma[4]-(tracegammam1/3.))*(Thetensorgamma[4]-(tracegammam1/3.))*Thetensoreta[4]+(Thetensorgamma[8]-(tracegammam1/3.))*(Thetensorgamma[8]-(tracegammam1/3.))*Thetensoreta[8])+(4.*constduekT)*((Thetensorgamma[1]*Thetensorgamma[1]/(invtraceeta+3.*Thetensoreta[8]))+(Thetensorgamma[2]*Thetensorgamma[2]/(invtraceeta+3.*Thetensoreta[4]))+(Thetensorgamma[5]*Thetensorgamma[5]/(invtraceeta+3.*Thetensoreta[0])));
                                  
                                     
    const double gdtr=g/30.;

  const double forcesquare=vecforce[0]*vecforce[0]+vecforce[1]*vecforce[1]+vecforce[2]*vecforce[2];
  for(int i=0; i<_dim; ++i){ 
    for(int j=0; j<_dim; ++j){
     double outd=0.;
     double diag=0.;
      outd=gdtr*vecforce[i]*vecforce[j];
      if(i==j){  
      	diag=(g/10.)*forcesquare+(k_B*constT*tracegammam1/3.);
	}
	veckeff.push_back(outd+diag);
	cout<<outd+diag<<" ";
         outfilew7<<outd+diag<<" ";
	 
    }
    cout<<endl;
       outfilew7<<endl;
 
  }
  return veckeff;
};

unsigned int Eulero::_GetRandom(){
	unsigned int ris;
	FILE *devran = fopen("/dev/urandom","r");
	fread(&ris, 4, 1, devran);
	fclose(devran);
	return ris;
}

vector<double> Eulero::_Random_Numb_Gen(){
    vector<double> _vec;
     const unsigned long long int myi=MAX_SIZE*3;
    if (myi>_vec.max_size()) {
        cout<<"MAX SIZE TOO BIG"<<endl;
        
        exit(-1);
    }
    _vec.resize(myi);
    unsigned int seed;
    
    seed=_GetRandom();
    default_random_engine generator (seed);
    normal_distribution<double> distribution(0.,1.0);
    for (unsigned long long int i=0; i<myi; ++i) {
        _vec[i]=distribution(generator);
    }
    
    return _vec;
};



void Eulero::Mean(vector<double> &xvec,vector<double> &yvec,vector<double> &zvec){
  
  


   double meanx, meany, meanz;
   double meanxsquare, meanysquare, meanzsquare;
   double meanxysquare, meanyzsquare, meanzxsquare;

   double meanxxx, meanyyy, meanzzz;
   double meanxyx, meanxyz, meanxzx;
   double meanyxy, meanyxz, meanyzy;
   double meanzxz, meanzxy, meanzzy;
 double meanxxxx, meanyyyy, meanzzzz;
        _vecmeanpassix.resize(PASSI );
    _vecmeanpassiy.resize(PASSI );
    _vecmeanpassiz.resize(PASSI );
    _vecmeanpassixsquare.resize(PASSI );
    _vecmeanpassiysquare.resize(PASSI );
    _vecmeanpassizsquare.resize(PASSI );
    _vecmeanpassixysquare.resize(PASSI );
    _vecmeanpassiyzsquare.resize(PASSI );
    _vecmeanpassizxsquare.resize(PASSI );
   _vecmeanpassixxx.resize(PASSI );
   _vecmeanpassiyyy.resize(PASSI );
   _vecmeanpassizzz.resize(PASSI );
   _vecmeanpassixyx.resize(PASSI );
   _vecmeanpassixyz.resize(PASSI );
   _vecmeanpassixzx.resize(PASSI );
   _vecmeanpassizzy.resize(PASSI );
   _vecmeanpassiyxy.resize(PASSI );
   _vecmeanpassiyxz.resize(PASSI );
   _vecmeanpassiyzy.resize(PASSI );
   _vecmeanpassizxz.resize(PASSI );
   _vecmeanpassizxy.resize(PASSI );
   _vecmeanpassixxxx.resize(PASSI );
   _vecmeanpassiyyyy.resize(PASSI );
   _vecmeanpassizzzz.resize(PASSI );
   
   unsigned  int i1=0;
   unsigned  int i2=0;
   unsigned  int i3=0;
  
   for(unsigned int j=0; j<PASSI; ++j){
   meanx=0.;
   meany=0.;
   meanz=0.;
   meanxsquare=0.;
   meanysquare=0.;
   meanzsquare=0.;
   meanxysquare=0.;
   meanyzsquare=0.;
   meanzxsquare=0.;
   meanxxx=0.;
   meanyyy=0.;
   meanzzz=0.;

   meanxxxx=0.;
   meanyyyy=0.;
   meanzzzz=0.;
   
   meanxyx=0.;
   meanxyz=0.;
   meanxzx=0.;
   meanyxy=0.;
   meanyxz=0.;
   meanyzy=0.;
   meanzxz=0.;
   meanzxy=0.;
   meanzzy=0.;
    for (unsigned int i=0; i<CONF; ++i) {
        i1=i*PASSI;
        meanx+=xvec[i1+j];//-xvec[i1+1];
        meany+=yvec[i1+j];//-yvec[i1+1];
        meanz+=zvec[i1+j];//-zvec[i1+1] ;
    }
       
    _vecmeanpassix[j] = (double)meanx/(double)CONF;
    _vecmeanpassiy[j] = (double)meany/(double)CONF;
    _vecmeanpassiz[j] = (double)meanz/(double)CONF;

      for (unsigned int i=0; i<CONF; ++i) {
          i2=i*PASSI;
        meanxsquare+=(xvec[i2+j]-_vecmeanpassix[j])*(xvec[i2+j]-_vecmeanpassix[j]);
        meanysquare+=(yvec[i2+j]-_vecmeanpassiy[j])*(yvec[i2+j]-_vecmeanpassiy[j]);
        meanzsquare+=(zvec[i2+j]-_vecmeanpassiz[j])*(zvec[i2+j]-_vecmeanpassiz[j]);
	meanxysquare+=(xvec[i2+j]-_vecmeanpassix[j])*(yvec[i2+j]-_vecmeanpassiy[j]);
	meanyzsquare+=(zvec[i2+j]-_vecmeanpassiz[j])*(yvec[i2+j]-_vecmeanpassiy[j]);
	meanzxsquare+=(xvec[i2+j]-_vecmeanpassix[j])*(zvec[i2+j]-_vecmeanpassiz[j]);
    }
      
   _vecmeanpassixsquare[j] = (double)meanxsquare/(double)(CONF);
   _vecmeanpassiysquare[j] = (double)meanysquare/(double)(CONF);
   _vecmeanpassizsquare[j] = (double)meanzsquare/(double)(CONF);
   _vecmeanpassixysquare[j]=(double)meanxysquare/(double)(CONF);
   _vecmeanpassiyzsquare[j]=(double)meanyzsquare/(double)(CONF);
   _vecmeanpassizxsquare[j]=(double)meanzxsquare/(double)(CONF);


   for (unsigned int i=0; i<CONF; ++i) { //here I should compute the skewness of the distribution.
          i3=i*PASSI;
        meanxxx+=(xvec[i3+j]-_vecmeanpassix[j])*(xvec[i3+j]-_vecmeanpassix[j])*(xvec[i3+j]-_vecmeanpassix[j]);
        meanyyy+=(yvec[i3+j]-_vecmeanpassiy[j])*(yvec[i3+j]-_vecmeanpassiy[j])*(yvec[i3+j]-_vecmeanpassiy[j]);
        meanzzz+=(zvec[i3+j]-_vecmeanpassiz[j])*(zvec[i3+j]-_vecmeanpassiz[j])*(zvec[i3+j]-_vecmeanpassiz[j]);
	meanxyx+=(xvec[i3+j]-_vecmeanpassix[j])*(yvec[i3+j]-_vecmeanpassiy[j])*(xvec[i3+j]-_vecmeanpassix[j]);
	meanxyz+=(zvec[i3+j]-_vecmeanpassiz[j])*(yvec[i3+j]-_vecmeanpassiy[j])*(xvec[i3+j]-_vecmeanpassix[j]);
	meanxzx+=(xvec[i3+j]-_vecmeanpassix[j])*(zvec[i3+j]-_vecmeanpassiz[j])*(xvec[i3+j]-_vecmeanpassix[j]);
        meanzzy+=(zvec[i3+j]-_vecmeanpassiz[j])*(yvec[i3+j]-_vecmeanpassiy[j])*(zvec[i3+j]-_vecmeanpassiz[j]);
	meanyxy+=(yvec[i3+j]-_vecmeanpassiy[j])*(yvec[i3+j]-_vecmeanpassiy[j])*(xvec[i3+j]-_vecmeanpassix[j]);
	meanyxz+=(xvec[i3+j]-_vecmeanpassix[j])*(zvec[i3+j]-_vecmeanpassiz[j])*(yvec[i3+j]-_vecmeanpassiy[j]);
	meanyzy+=(zvec[i3+j]-_vecmeanpassiz[j])*(yvec[i3+j]-_vecmeanpassiy[j])*(yvec[i3+j]-_vecmeanpassiy[j]);
	meanzxz+=(zvec[i3+j]-_vecmeanpassiz[j])*(zvec[i3+j]-_vecmeanpassiz[j])*(xvec[i3+j]-_vecmeanpassix[j]);
	meanzxy+=(xvec[i3+j]-_vecmeanpassix[j])*(zvec[i3+j]-_vecmeanpassiz[j])*(yvec[i3+j]-_vecmeanpassiy[j]);
    }

   _vecmeanpassixxx[j] = (double)meanxxx/(double)(CONF);
   _vecmeanpassiyyy[j] = (double)meanyyy/(double)(CONF);
   _vecmeanpassizzz[j] = (double)meanzzz/(double)(CONF);
   _vecmeanpassixyx[j] = (double)meanxyx/(double)(CONF);
   _vecmeanpassixyz[j] = (double)meanxyz/(double)(CONF);
   _vecmeanpassixzx[j] = (double)meanxzx/(double)(CONF);
   _vecmeanpassizzy[j] = (double)meanzzy/(double)(CONF);
   _vecmeanpassiyxy[j] = (double)meanyxy/(double)(CONF);
   _vecmeanpassiyxz[j] = (double)meanyxz/(double)(CONF);
   _vecmeanpassiyzy[j] = (double)meanyzy/(double)(CONF);
   _vecmeanpassizxz[j] = (double)meanzxz/(double)(CONF);
   _vecmeanpassizxy[j] = (double)meanzxy/(double)(CONF);

   
    for (unsigned int i=0; i<CONF; ++i) { //here I should compute the skewness of the distribution.
          i3=i*PASSI;
        meanxxxx+=(xvec[i3+j]-_vecmeanpassix[j])*(xvec[i3+j]-_vecmeanpassix[j])*(xvec[i3+j]-_vecmeanpassix[j])*(xvec[i3+j]-_vecmeanpassix[j]);
        meanyyyy+=(yvec[i3+j]-_vecmeanpassiy[j])*(yvec[i3+j]-_vecmeanpassiy[j])*(yvec[i3+j]-_vecmeanpassiy[j])*(yvec[i3+j]-_vecmeanpassiy[j]);
        meanzzzz+=(zvec[i3+j]-_vecmeanpassiz[j])*(zvec[i3+j]-_vecmeanpassiz[j])*(zvec[i3+j]-_vecmeanpassiz[j])*(zvec[i3+j]-_vecmeanpassiz[j]);

    }

    _vecmeanpassixxxx[j] = (double)meanxxxx/(double)(CONF);
    _vecmeanpassiyyyy[j] = (double)meanyyyy/(double)(CONF);
    _vecmeanpassizzzz[j] = (double)meanzzzz/(double)(CONF);
   }
   
};

void Eulero::Integration(vector<double> &posizionextempot, vector<double> &posizioneytempot, vector<double> &posizioneztempot, Tensor &Thetensorgamma, Tensor &Thetensoreta, Temperature &TheTemperature, Force &TheForce, Torque &TheTorque){

    // Check this point
    _mystring="transandrotinLabFrame";  
 
    Quaternion TheQuaternion;
    Quaternion TheQuaternionc;
    Quaternion TheQuaternion1;
    Quaternion TheQuaternion1c;
    Quaternion TheQuaternionOmega;
    Quaternion TheQuaterniontrans;
    Quaternion TheQuaterniontransm1;
    Quaternion TheQuaternionForce;
    Quaternion TheQuaternionTorqueLab;

    double zero=0.;
    _constT=TheTemperature.ConstTemperature();
    TheTorque.ConstTorque();
    TheForce.ConstForce();  
    /* Rotational vector*/
    vector<double> vectemprot;
    vector<double> vectemprotsquare;
    vector<double> _vectorpositionrot;
    vector<double> vecforcelab;
    vector<double> vectorquelab;
    vecforcelab=TheForce.ConstForce();
    vectorquelab=TheTorque.ConstTorque();
    
    TheQuaternionForce.PushQuaternion(zero, vecforcelab);
    TheQuaternionTorqueLab.PushQuaternion(zero, vectorquelab);

    _vectorpositionrot.resize(3);
    vectemprot.resize(3);
    vectemprotsquare.resize(3);
    
    _mystringgammax<<Thetensorgamma.ElementTensor(0, 0);
    _mystringgammay<<Thetensorgamma.ElementTensor(1, 1);
    _mystringgammaz<<Thetensorgamma.ElementTensor(2, 2);
    _mystringgammaxy<<Thetensorgamma.ElementTensor(0, 1);
    _mystringgammaxz<<Thetensorgamma.ElementTensor(0, 2);
    _mystringgammayz<<Thetensorgamma.ElementTensor(1, 2);
    _mystringforcex<<vecforcelab[0];
    _mystringforcey<<vecforcelab[1];
    _mystringforcez<<vecforcelab[2];

    KeffRS(Thetensorgamma,Thetensoreta,TheTemperature,TheForce);
    vector<double> vectemp;
    vector<double> vectemppossquare;
    const double invercetraceeta=Thetensoreta.Trace();
    double x=0., y=0. ,z=0.;
    vectemp.resize(3);
    Quaternion q1;
    Quaternion q5;
    Quaternion q6;
    Quaternion q2;
    Quaternion q3;
    Quaternion q4;
    vector<double> vectorquebody;
    vector<double> vecforcebody;
    vector<double> vectemp1;
    vector<double> sqrtThetensorgamma;
    vectemp1.resize(3);
    vectorquebody.resize(3);
    vecforcebody.resize(3);
    sqrtThetensorgamma=Thetensorgamma.InvVec();
   /* unsigned int seed, seed1, seed2, seedr, seedr1, seedr2;
    
    seed=_GetRandom();
    seed1=_GetRandom();
    seed2=_GetRandom();
    seedr=_GetRandom();
    seedr1=_GetRandom();
    seedr2=_GetRandom();
    default_random_engine generator (seed);
    default_random_engine generator1 (seed1);
    default_random_engine generator2 (seed2);
    default_random_engine generatorrot (seedr);
    default_random_engine generatorrot1 (seedr1);
    default_random_engine generatorrot2 (seedr2);
    normal_distribution<double> distribution(0.,1.0);
    normal_distribution<double> distribution1(0.,1.0);
    normal_distribution<double> distribution2(0.,1.0);
    normal_distribution<double> distributionrot(0.,1.0);
    normal_distribution<double> distributionrot1(0.,1.0);
    normal_distribution<double> distributionrot2(0.,1.0);*/
    vector<double> vecnoisetrans;
    vector<double> vecnoiserot;
    vector<double> vecnoisetrans1;
    vector<double> vecnoiserot1;
    vector<double> vecnoisetrans2;
    vector<double> vecnoiserot2;
    double mynewconst;
    unsigned int mycount, mycount1;
    mynewconst=-k_B*_constT*dt*invercetraceeta;
    for (unsigned int conf=0; conf<CONF; ++conf) {
           cout<<conf<<endl;
       TheQuaternion.RotationQuaternion();
        vectemp[0]=0.;
        vectemp[1]=0.;
        vectemp[2]=0.;
        vecnoisetrans=_Random_Numb_Gen();
        vecnoiserot=_Random_Numb_Gen();
        vecnoisetrans1=_Random_Numb_Gen();
        vecnoiserot1=_Random_Numb_Gen();
        vecnoisetrans2=_Random_Numb_Gen();
        vecnoiserot2=_Random_Numb_Gen();
	mycount=(deltat/dt);
	mycount1=1;
       for (unsigned int time=1; time<MAX_SIZE; ++time) {
        TheQuaternionc=TheQuaternion.conjugate();
        q5=TheQuaternionc*TheQuaternionTorqueLab;
        q6=q5*TheQuaternion;
        vectorquebody[0]=q6.xQuaternion();
        vectorquebody[1]=q6.yQuaternion();
        vectorquebody[2]=q6.yQuaternion();
           
         for (unsigned int j=0; j<_dim; ++j) {
             
             _vectorpositionrot[0]+=(Thetensoreta[j])*vectorquebody[j]*dt+sqrt(2.*k_B*_constT*dt)*(sqrt(Thetensoreta[j]))*vecnoiserot[3*time+j];
             _vectorpositionrot[1]+=(Thetensoreta[_dim+j])*vectorquebody[j]*dt+sqrt(2.*k_B*_constT*dt)*(sqrt(Thetensoreta[_dim+j]))*vecnoiserot1[3*time+j];
             _vectorpositionrot[2]+=(Thetensoreta[2*_dim+j])*vectorquebody[j]*dt+sqrt(2.*k_B*_constT*dt)*(sqrt(Thetensoreta[2*_dim+j]))*vecnoiserot2[3*time+j];
                 
            }
		  TheQuaternionOmega.PushQuaternion(mynewconst, _vectorpositionrot);
		  TheQuaternion1=TheQuaternion*TheQuaternionOmega;
		  TheQuaternion1.Prodforonehalf(0.5);
		  TheQuaternion1c=TheQuaternion1.conjugate();
		 
           
		  _vectorpositionrot[0]=0.;
		  _vectorpositionrot[1]=0.;
		  _vectorpositionrot[2]=0.;
	 
	   
        q2=TheQuaternionc*TheQuaternionForce;
        q3=q2*TheQuaternion;
        vecforcebody[0]=q3.xQuaternion();
        vecforcebody[1]=q3.yQuaternion();
        vecforcebody[2]=q3.zQuaternion();
	   
                for (int j=0; j<_dim; ++j) {

		            _vectorposition[0]+=(Thetensorgamma[j])*vecforcebody[j]*dt+sqrt(2.*k_B*_constT*dt)*(sqrtThetensorgamma[j])*vecnoisetrans[3*time+j];
                    _vectorposition[1]+=(Thetensorgamma[_dim+j])*vecforcebody[j]*dt+sqrt(2.*k_B*_constT*dt)*(sqrtThetensorgamma[_dim+j])*vecnoisetrans1[3*time+j];
                    _vectorposition[2]+=(Thetensorgamma[2*_dim+j])*vecforcebody[j]*dt+sqrt(2.*k_B*_constT*dt)*(sqrtThetensorgamma[2*_dim+j])*vecnoisetrans2[3*time+j];
	
		}
           vectemp=_vectorposition;
           TheQuaterniontrans.PushQuaternion(zero, vectemp);
           
           q1=TheQuaternion*TheQuaterniontrans;
           q4=q1*TheQuaternionc;
           vectemp1[0]= q4.xQuaternion();
           vectemp1[1]= q4.yQuaternion();
           vectemp1[2]= q4.zQuaternion();
           
           x+=vectemp1[0];
           y+=vectemp1[1];
           z+=vectemp1[2];
           
           TheQuaternion=TheQuaternion+TheQuaternion1;
           TheQuaternion.NormalizationQuaternion();
	   if( mycount==time ){
           posizionextempot[conf*PASSI+mycount1]=posizionextempot[conf*PASSI+mycount1-1]+x;
           posizioneytempot[conf*PASSI+mycount1]=posizioneytempot[conf*PASSI+mycount1-1]+y;
           posizioneztempot[conf*PASSI+mycount1]=posizioneztempot[conf*PASSI+mycount1-1]+z;
           ++mycount1;
           mycount=mycount1*(deltat/dt);
           x=0.;
           y=0.;
           z=0.;
	   }
            _vectorposition[0]=0.;
            _vectorposition[1]=0.;
            _vectorposition[2]=0.;
           if(mycount1==the_time)break;
        
       }
    
       }
     Mean(posizionextempot,posizioneytempot,posizioneztempot);
     PrintonFileTrajectory(posizionextempot,posizioneytempot,posizioneztempot);
 


};

void Eulero::PrintonFile(){
    string mystring;
    ostringstream T_str;
    T_str<<_constT;
    const string title="Possquare";
    const string dat=".txt";
    mystring=title;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();
    mystring+=T_str.str();
    mystring+=dat;
    ofstream outfilew(mystring.c_str(), ios_base::app );
    for (int time=1; time<PASSI; ++time) {
        outfilew<<time*deltat<<" "<<_vecmeanpassixsquare[time]<<" "<<_vecmeanpassiysquare[time]<<" "<<_vecmeanpassizsquare[time]<<endl;
    }
    // Please you should check the print function, I do not know if everything works well 
    
    string mystring1;
    const string title1="Diffusion";
    const string dat1=".txt";
    mystring1=title1;
    mystring1+=_mystring;
    mystring1+=_mystringgammax.str();
    mystring1+=_mystringgammay.str();
    mystring1+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring1+=_mystringforcex.str();
    mystring1+=_mystringforcey.str();
    mystring1+=_mystringforcez.str();
    mystring1+=T_str.str();
    mystring1+=dat1;
    ofstream outfilew1(mystring1.c_str(), ios_base::app );
    
    for (int time=1; time<PASSI; ++time) {
      outfilew1<<time*deltat<<" "<<(0.5)*_vecmeanpassixsquare[time]/(double)((double)time*deltat)<<" "<<(0.5)*_vecmeanpassiysquare[time]/(double)((double)time*deltat)<<" "<<(0.5)*_vecmeanpassizsquare[time]/(double)((double)time*deltat)<<" "<<(0.5)*_vecmeanpassixysquare[time]/(double)((double)time*deltat)<<" "<<(0.5)*_vecmeanpassiyzsquare[time]/(double)((double)time*deltat)<<" "<<(0.5)*_vecmeanpassizxsquare[time]/(double)((double)time*deltat)<<endl;

    }
    
    const string title2="BallisticVelocity";
    const string dat2=".txt";
    mystring=title2;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();    
    mystring+=T_str.str();
    mystring+=dat2;
    ofstream outfilew2(mystring.c_str(), ios_base::app );
    for (int time=1; time<PASSI; ++time) {
         outfilew2<<time*deltat<<" "<< _vecmeanpassix[time]/(double)((double)time*deltat)<<" "<< _vecmeanpassiy[time]/(double)((double)time*deltat)<<" "<< _vecmeanpassiz[time]/(double)((double)time*deltat)<<endl;
      // Da cambiare dobbiamo sicuramente diminuire la memoria
    }
    const string titleskew="3momentdivtdt";
    const string datskew=".txt";
    mystring=titleskew;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();    
    mystring+=T_str.str();
    mystring+=datskew;
    

    ofstream outfilewskew(mystring.c_str(), ios_base::app );
    for (int time=1; time<PASSI; ++time) {
      outfilewskew<<time*deltat<<" "<<  _vecmeanpassixxx[time]/(double)((double)time*deltat)<<" "<<_vecmeanpassiyyy[time]/(double)((double)time*deltat)<<" "<<_vecmeanpassizzz[time]/(double)((double)time*deltat)<<" "<<  _vecmeanpassixyx[time]/(double)((double)time*deltat)<<" "<<  _vecmeanpassixyz[time]/(double)((double)time*deltat)<<" "<<  _vecmeanpassixzx[time]/(double)((double)time*deltat)<<" "<<  _vecmeanpassizzy[time]/(double)((double)time*deltat)<<" "<<  _vecmeanpassiyxy[time]/(double)((double)time*deltat)<<" "<<  _vecmeanpassiyxz[time]/(double)((double)time*deltat)<<" "<<  _vecmeanpassiyzy[time]/(double)((double)time*deltat)<<" "<<  _vecmeanpassizxz[time]/(double)((double)time*deltat)<<" "<<  _vecmeanpassizxy[time]/(double)((double)time*deltat)<<endl;
      // Da cambiare dobbiamo sicuramente diminuire la memoria
    }

      const string titleskew1="3moment";
    const string datskew1=".txt";
    mystring=titleskew1;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();    
    mystring+=T_str.str();
    mystring+=datskew1;
    

    ofstream outfilewskew1(mystring.c_str(), ios_base::app );
    for (int time=1; time<PASSI; ++time) {
      outfilewskew1<<time*deltat<<" "<<  _vecmeanpassixxx[time]<<" "<<_vecmeanpassiyyy[time]<<" "<<_vecmeanpassizzz[time]<<" "<<  _vecmeanpassixyx[time]<<" "<<  _vecmeanpassixyz[time]<<" "<<  _vecmeanpassixzx[time]<<" "<<  _vecmeanpassizzy[time]<<" "<<  _vecmeanpassiyxy[time]<<" "<<  _vecmeanpassiyxz[time]<<" "<<  _vecmeanpassiyzy[time]<<" "<<  _vecmeanpassizxz[time]<<" "<<  _vecmeanpassizxy[time]<<endl;
      // Da cambiare dobbiamo sicuramente diminuire la memoria
    }



      const string title4mom="4moment";
    const string dat4mom=".txt";
    mystring=title4mom;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();    
    mystring+=T_str.str();
    mystring+= dat4mom;
    

    ofstream outfile4mom(mystring.c_str(), ios_base::app );
    for (int time=1; time<PASSI; ++time) {
      outfile4mom<<time*deltat<<" "<<  _vecmeanpassixxxx[time]<<" "<<_vecmeanpassiyyyy[time]<<" "<<_vecmeanpassizzzz[time]<<" "<<  _vecmeanpassixxxx[time]-(3*(_vecmeanpassixsquare[time])*(_vecmeanpassixsquare[time]))<<" "<<_vecmeanpassiyyyy[time]-(3*(_vecmeanpassiysquare[time])*(_vecmeanpassiysquare[time]))<<" "<<_vecmeanpassizzzz[time]-(3*(_vecmeanpassizsquare[time])*(_vecmeanpassizsquare[time]))<<endl;
      // Da cambiare dobbiamo sicuramente diminuire la memoria
    }

         const string titleskew2="skewness";
    const string datskew2=".txt";
    mystring=titleskew2;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();    
    mystring+=T_str.str();
    mystring+=datskew2;
    

    ofstream outfilewskew2(mystring.c_str(), ios_base::app );
    for (int time=1; time<PASSI; ++time) {
      outfilewskew2<<time*deltat<<" "<<  _vecmeanpassixxx[time]/(_vecmeanpassixsquare[time]*sqrt(_vecmeanpassixsquare[time]))<<" "<<_vecmeanpassiyyy[time]/(_vecmeanpassiysquare[time]*sqrt(_vecmeanpassiysquare[time]))<<" "<<_vecmeanpassizzz[time]/(_vecmeanpassizsquare[time]*sqrt(_vecmeanpassizsquare[time]))<<endl;
      // Da cambiare dobbiamo sicuramente diminuire la memoria
    }
    
};

void Eulero::PrintonFileTrajectory(vector<double> & xvec,vector<double> & yvec,vector<double> & zvec){
    string mystring;
    ostringstream T_str;
    T_str<<_constT;
    const string title2="Trajectories";
    const string dat2=".txt";
    mystring=title2;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();
    mystring+=T_str.str();
    mystring+=dat2;
    ofstream outfilew3(mystring.c_str(), ios_base::app );
    for(int i=0; i<1; ++i){
    for (int time=0; time<PASSI; ++time) {
         outfilew3<<time*deltat<<" "<< xvec[i*PASSI+time]<<" "<< yvec[i*PASSI+time]<<" "<< zvec[i*PASSI+time]<<endl;
      // Da cambiare dobbiamo sicuramente diminuire la memoria
    }
    }
    


    const string title3="Trajectories2";
    const string dat3=".txt";
    mystring=title3;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();
    mystring+=T_str.str();
    mystring+=dat3;
    ofstream outfilew4(mystring.c_str(), ios_base::app );
    for(int i=1; i<2; ++i){
        for (int time=0; time<PASSI; ++time) {
            outfilew4<<time*deltat<<" "<< xvec[i*PASSI+time]<<" "<< yvec[i*PASSI+time]<<" "<< zvec[i*PASSI+time]<<endl;
         
        }
    }
    
    const string title4="Trajectories3";
    const string dat4=".txt";
    mystring=title4;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();
    mystring+=T_str.str();
    mystring+=dat4;
    ofstream outfilew5(mystring.c_str(), ios_base::app );
    for(int i=2; i<3; ++i){
        for (int time=0; time<PASSI; ++time) {
            outfilew5<<time*deltat<<" "<< xvec[i*PASSI+time]<<" "<< yvec[i*PASSI+time]<<" "<< zvec[i*PASSI+time]<<endl;
            
        }
    }

    const string title5="pdfasymmetry";
    const string dat5=".txt";
    mystring=title5;
    mystring+=_mystring;
    mystring+=_mystringgammax.str();
    mystring+=_mystringgammay.str();
    mystring+=_mystringgammaz.str();
    mystring+=_mystringgammaxy.str();
    mystring+=_mystringgammayz.str();
    mystring+=_mystringgammaxz.str();
    mystring+=_mystringforcex.str();
    mystring+=_mystringforcey.str();
    mystring+=_mystringforcez.str();
    mystring+=T_str.str();
    mystring+=dat5;
    ofstream outfilew8(mystring.c_str(), ios_base::app );
    for(int i=0; i<CONF; ++i){
       
      outfilew8<< xvec[i*PASSI+(PASSI-1)]<<" "<< yvec[i*PASSI+(PASSI-1)]<<" "<< zvec[i*PASSI+(PASSI-1)]<<endl;
            
        
    }

};

#endif
