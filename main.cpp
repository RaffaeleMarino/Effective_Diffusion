//
//  main.cpp
//  FinalEffDiff
//
//  Created by Raffaele Marino on 12/12/15.
//  Copyright (c) 2015 Raffaele Marino. All rights reserved.
//

#include "header_s.h"
#include "Class_Quaternion.hpp"
#include "Class_Tensor.hpp"
#include "Class_Temperature.hpp"
#include "Class_EULER_RK_INTEGR.hpp"
#include "Class_Force.hpp"
#include "Class_Torque.hpp"




void help();

int main(int argc, char * argv[])
{
    Tensor TheEta;
    Tensor TheGamma;
    Force TheForce;
    Torque TheTorque;
    Temperature TheTemperature;
    int c;
    int dim;
	while ((c = getopt(argc, argv, "w:h?")) != -1)
		switch (c) {
			case 'w':
                dim=atoi(argv[2]);
                switch (dim) {
                    case 2:
                        /*The implementation of the code for 2-d will be done soon*/
                        cout<<"The code has to be implemented for 2-dimension."<<endl;
                        cout<<"You can use this 3-d with a coordinate equal to zero"<<endl;
                        exit(-1);
                        
                        
                        break;
                        
                    case 3:
                        /*Important here the tensor is the inverse tensor*/
                        cout<<"You have chosen 3-dimension"<<endl;
                        cout<<"INVERSE Translational Friction Tensor"<<endl;
                        TheGamma.NewTensorBodysTRANS(argv, dim);
                        cout<<"INVERSE Rotational Friction Tensor"<<endl;
                        TheEta.NewTensorBodysROT(argv, dim);
                        TheForce.Inizialization(argv, dim);
                        TheTorque.Inizialization(argv, dim);
                        TheTemperature.Inizialization(argv, dim);
                        break;
                        
                    default:
                        help();
                        exit(-1);
                        break;
                }
                
				break;
				
			case 'h':
			default:
				fprintf(stderr, "%s [options]\n"
						"\n"
						"-w= input values for Tensors, Force, Torque and Temperature \n"
						,argv[0]);
				help();
				exit(-1);
				break;
		}
    

	Eulero TheEulerTransandRot(dim);
    vector<double> posizionextempot;
	vector<double> posizioneytempot;
	vector<double> posizioneztempot;
	posizionextempot.resize(CONF*PASSI);
	posizioneytempot.resize(CONF*PASSI);
	posizioneztempot.resize(CONF*PASSI);
	TheEulerTransandRot.Integration(posizionextempot, posizioneytempot, posizioneztempot, TheGamma,TheEta, TheTemperature, TheForce, TheTorque);
	TheEulerTransandRot.PrintonFile();
    return 0;
}

void help(){
    cout << ""<<endl;
    cout << ""<<endl;
    cout << ""<<endl;
    cout << ""<<endl;
	cout << "***********************************************This is the help function***********************************************"<<endl;
    cout << ""<<endl;
    cout << ""<<endl;
    cout << "********************************************************************************************************************************************************************************************"<<endl;
    cout << ""<<endl;
    cout << ""<<endl;
    cout << "This code was written by Raffaele Marino. If you need information you can contact the author: raffaele.marino@nordita.org."<<endl;
    cout << ""<<endl;
    cout << ""<<endl;
    cout << "********************************************************************************************************************************************************************************************"<<endl;
    cout << ""<<endl;
    cout << ""<<endl;
    cout << "If you want to use this program you should follow the istructions written below:"<<endl;
    cout << ""<<endl;
    cout << "if you use this program you should know that:"<<endl;
    cout << ""<<endl;
    cout << "it can work in 2-dimension and 3-dimension."<<endl;
    cout << ""<<endl;
    cout <<"All the tensors are symmeric."<<endl;
    cout << ""<<endl;
    cout << "Translational Friction Tensor is written in Body's System. Diagonal Values of the Tensor have to be set to zero"<<endl;
    cout << ""<<endl;
    cout << "Rotational Friction Tensor is written in Body's System."<<endl;
    cout << ""<<endl;
    cout << "Please, if you want to work in 2-dimension you should not write the elements of tensors and the components of the vectors out of range."<<endl;
    cout << ""<<endl;
    cout << "You shold write the components of the Force and Torque vectors as a constant."<<endl;
    cout << ""<<endl;
    cout << "You shold write the Temperature as a constant."<<endl;
    cout << ""<<endl;
    cout <<"If you want to use formulae for the Torque, Force or Temperature, you should contact the author: raffaele.marino@nordita.org."<<endl;
    cout << ""<<endl;
    cout << ""<<endl;
    cout << "********************************************************************************************************************************************************************************************"<<endl;
    cout << ""<<endl;
    cout << ""<<endl; // It is important change this part of the code.
    cout << "If you use -w: [-w ] [-dimension] [Diagonal Value of Translational Friction Tensor (1,1)] [Diagonal Value of Translational Friction Tensor (2,2)] [Diagonal Value of Translational Friction Tensor (3,3)] [OUT Diagonal Value of Translational Friction Tensor (1,2)] [OUT Diagonal Value of Translational Friction Tensor (1,3)] [OUT Diagonal Value of Translational Friction Tensor (2,3)][Diagonal Values of Rotational Friction Tensor (1,1)] [Diagonal Values of Rotational Friction Tensor (2,2)] [Diagonal Values of Rotational Friction Tensor (3,3)] [Force component vector 1] [Force component vector 2][Force component vector 3] [Torque component vector 1][Torque component vector 2][Torque component vector 3] [Temperature]"<<endl;
}


