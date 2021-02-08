#ifndef APControl_H
#define APControl_H

//#define DEBUG_MODE
// Valores por defecto
#define MAX_INC_CTRL  0
//Límites absolutos de la acción de control
#define MAX_CTRL 1.0    
#define MIN_CTRL 0
#define DELTA_b 0 // Con delta_b = 0, caso ideal
#define ADAPT_FACT 0.2

#include<Arduino.h>

/*
*		La clase transferFunction define una funcion de transferencia discreta que viene dada 
		por los parámetros de la ecuación en diferencias correspondiente:
		y(k) = a1·y(k-1) + a2·y(k-2) +...+ an·y(k-n) + b1·u(k-1) + b2·u(k-2) + ... + bn·u(k-n)
		Donde a1, a2, ...,an son los parámetros de salida y u1, u2, ..., un son los parámetros
		de entrada
*/
class transferFunction {
	public:
		transferFunction(float *o, byte o_l, float *i, byte i_l);
	
		float *output;
		float *input;
		byte outputLen;
		byte inputLen;
};
/*
*		La clase APControl implementa un controlador Adaptativo Predictivo
*/
class APControl{
  


//Parametros del modelo Adaptativo Predictivo
	transferFunction *model;

//Parámetros de la trayectoria deseada 
	transferFunction *traject;

//Variables del proceso

	float *y;		//La variable de salida
	float *u;		//La variable de entrada
	float *ySP;		// El set-point del sistema
	float yDNext;   // La salida deseada en k+1
	float err;      // El error de estimación a posteriori
	byte lambda;		//Horizonte de prediccion
	float maxControl; //Máxima acción de control
	float minControl; //Mínima acción de control
	float maxIncControl; //Máximo incremento absoluto de control
	float deltaB;	
	float adaptFact;

//Coeficientes del control predictivo extendido
	float **e;
	float **g;
	float h;
	float **fi;
	float **delta;
	float mi;
	float *ni;
	float *gamma;

	void extendedCoefs();
	
  public:

	APControl( transferFunction *mod, transferFunction *traj, byte l); //Estrategia extendida
	APControl( transferFunction *mod, transferFunction *traj); //Estrategia básica

	void setSetPoint(float sp);
	void setLambda(byte l);
    void setProcessOutput(float out);
    void setControlAction(float action);
	void setMaxControl(float m);
	void setMinControl(float m);
	void setMaxIncControl(float m);
	void setDeltaB(float db);
	void setAdaptFact(float adapt);
	void setModel(transferFunction *model);
	void setTraject(transferFunction *traject);


    float getYd();
    float getControlAction();
    float getPostErr();
	int getLambda();
    float getProcessOutput();
    float getSetPoint();
    void updateVars();
    void predictiveAction();
	void basicPredictiveAction();
    void adaptativeMechanism();


	float getAm1();
	float getAm2();
	float getAm3();
	float getBm1();
	float getBm2();
	float getBm3();
};

#endif
