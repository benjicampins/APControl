//////////////////////////////////////////////////////////
// CONTROL ADAPTATIVO PREDICTIVO PARA ARDUINO           //
// Estrategia extendida, formulación incremental
// Autor: Benjamín Campins                              //
// Versión: 7                                         //
// Fecha: 28/12/2020                                    //
//                                                      //
/////////////////////////////////////////////////////////

#include"Arduino.h"
#include"APControl.h"



/*
*				Función de transferencia de sistema discreto
				Definida por los parámetros de entrada b1, b2,...,bm
				y los parámetros de salida a1, a2, an. Siendo la ecuación en diferencias:
				y(k) = a1·y(k-1) + a2·y(k-2) + ... + an·y(k-n) + b1·u(k-1) + b2·u(k-2) + ... + bm·u(k-m)
				Donde "u" es la entrada, e "y" es la salida.
*/

transferFunction::transferFunction(float *o, byte oL, float *i, byte iL){
	output = o;
	input = i;
	outputLen = oL;
	inputLen = iL;
}
/*
*				Se deben pasar punteros a dos funciones de transferencia definidas por 
				el usuario, model representa el modelo adaptativo-predictivo, y
				traject, la trayectoria de referencia.
				El parámetro lambda es el horizonte de predicción, lamda=1 para la estrategia básica.
*/

APControl::APControl(transferFunction *mod, transferFunction *traj, byte lam){
	
	maxControl = MAX_CTRL;
	minControl = MIN_CTRL;
	maxIncControl = MAX_INC_CTRL;
	deltaB = DELTA_b;
	adaptFact = ADAPT_FACT;
	model = mod;
	traject = traj;
	y = new float[(mod->outputLen) + 1];
	u = new float[(mod->inputLen) + 1];
	ySP = new float[traj->inputLen];
	ni = new float[mod->outputLen];
	gamma = new float[mod->inputLen];
	lambda = lam;
	//Crea las matrices E, G, lambda y fi e introduce ceros en los valores
	e = new float *[model->outputLen+1];
    for(int i = 0; i<=(model->outputLen)+1;i++){
      e[i] = new float[lambda]{};
    }
	g = new float *[model->inputLen+1];
    for(int i = 0; i<=(model->inputLen)+1;i++){
      g[i] = new float[lambda]{};
    }
	fi = new float *[traject->outputLen+1];
    for(int i = 0; i<=(traject->outputLen)+1;i++){
      fi[i] = new float[lambda]{};
    }
	delta = new float *[traject->inputLen+1];
    for(int i = 0; i<=(traject->inputLen)+1;i++){
      delta[i] = new float[lambda]{};
    }
}


APControl::APControl(transferFunction *mod, transferFunction *traj){
	
	maxControl = MAX_CTRL;
	minControl = MIN_CTRL;
	maxIncControl = MAX_INC_CTRL;
	deltaB = DELTA_b;
	adaptFact = ADAPT_FACT;
	model = mod;
	traject = traj;
	y = new float[(mod->outputLen) + 1];
	u = new float[(mod->inputLen) + 1];
	ySP = new float[traj->inputLen];
	lambda = 1;
}

/*
*		Actaliza los valores de las variables, para avanzar al siguiente periodo de muestreo
*/
void APControl::updateVars(){
  int i;
  for(i=model->outputLen; i>0; i--) y[i] = y[i-1];
  for(i=model->inputLen; i>0; i--) u[i] = u[i-1];
  for(i=traject->inputLen ; i>0; i--) ySP[i] = ySP[i-1];
  
}
/*
		Calcula la trayectoria de referencia y la acción de control predictivo

*		Estrategia básica, lambda = 1
*/

void APControl::basicPredictiveAction(){
	
	int i;
	//Calculo salida deseada en k+1	
	yDNext = 0;
	for(i=0;i<traject->outputLen;i++) yDNext += traject->output[i]*y[i];
	for(i=0;i<traject->inputLen;i++) yDNext += traject->input[i]*ySP[i];
	//Calculo de la acción de control teórica
	float incU = yDNext - y[0];
	for(i=0; i<model->outputLen; i++) incU -= (model->output[i]*(y[i] - y[i+1]));
	for(i=1; i<model->inputLen; i++) incU -= (model->input[i]*(u[i] - u[i+1]));
	incU  /= model->input[0];
	//Comprobación de que la acción no supere el máximo incremento permitido
	if(incU > maxIncControl && maxIncControl != 0  ) incU = maxIncControl;
	else if(incU < (-maxIncControl) && maxIncControl != 0) incU = -maxIncControl;
	//Comprobación de que la acción de control está dentro de los limites absolutos
	if((u[1]+incU) < minControl) u[0] = minControl;
	else if(u[1]+incU > maxControl) u[0] = maxControl;
	//Si cumple todas las restricciones, aplica el incremento de la acción calculada
	else u[0] = u[1] + incU;
	
}

/*
	Calcula la trayectoria de referencia y la acción de control predictivo
	estrategia extendida
*/

void APControl::predictiveAction(){
	
	int i;
	extendedCoefs();
	//Calculo trayectoria de referencia en k+lambda
	yDNext = 0;
	for(i=0;i<traject->outputLen;i++) yDNext += fi[i][lambda-1]*y[i];
	for(i=1;i<traject->inputLen;i++) yDNext += delta[i][lambda-1]*ySP[i];
	yDNext += mi*ySP[0];
	//Calculo de la acción de control teórica
	float incU = yDNext - y[0];
	for(i=0; i<model->outputLen; i++) incU -= (ni[i]*(y[i] - y[i+1]));
	for(i=1; i<model->inputLen; i++) incU -= (gamma[i]*(u[i] - u[i+1]));
	incU  /= h;
	//Comprobación de que la acción no supere el máximo incremento permitido
	if(incU > maxIncControl && maxIncControl != 0  ) incU = maxIncControl;
	else if(incU < (-maxIncControl) && maxIncControl != 0) incU = -maxIncControl;
	//Comprobación de que la acción de control está dentro de los limites absolutos
	if((u[1]+incU) < minControl) u[0] = minControl;
	else if(u[1]+incU > maxControl) u[0] = maxControl;
	//Si cumple todas las restricciones, aplica el incremento de la acción calculada
	else   u[0] =  u[1] + incU;

}

/*
*		Ajusta los parámetros del modelo Adaptativo Predictivo
*/
void APControl::adaptativeMechanism(){
  float fiQuad = 0;
  int i = 0;
  for(i=0;i<model->outputLen; i++) fiQuad += y[i+1]*y[i+1];
  for(i=0;i<model->inputLen; i++) fiQuad += u[i+1]*u[i+1];
  float amNext[model->outputLen]; // Próximos valores del modelo
  float bmNext[model->inputLen];
  float prError = y[0];
  //Ecuación 5.16;
  for(i=0;i<model->outputLen; i++) prError -= model->output[i]*y[i+1];
  for(i=0;i<model->inputLen; i++) prError -= model->input[i]*u[i+1];
  //Ecuación 5.20
  float postError = prError / (1+fiQuad);
  for(i=0;i<model->outputLen; i++) amNext[i] = 
  model->output[i] + postError*y[i+1]*adaptFact;
  for(i=0;i<model->inputLen; i++) bmNext[i] = 
  model->input[i] + postError*u[i+1]*adaptFact;
  // Condiciones 6.14 
  if (abs(prError) >= deltaB){
    for(i=0;i<model->outputLen; i++) model->output[i] = amNext[i];
    for(i=0;i<model->inputLen; i++) model->input[i] = bmNext[i];
	}
	#ifdef DEBUG_MODE
	else Serial.println("No aplica adaptación!");
	#endif
	
 err = abs(postError);
}

void APControl::extendedCoefs(){
	int i, j;
	h = 0; 	mi = 0;
	//Ecuaciones (4.9)
	for(i=0;i<model->outputLen;i++)	e[i][0] = model->output[i];
	for(i=0;i<model->inputLen;i++)	g[i][0] = model->input[i];	
	for(j=1;j<lambda;j++) e[model->outputLen][j-1] = 0;
	for(j=1;j<lambda;j++) g[model->inputLen][j-1] = 0;
	//Ecuaciones (4.8)
	for(j=1;j<lambda;j++){
		for(i=0;i<model->outputLen;i++){
			e[i][j] = e[0][j-1]*model->output[i] + e[i+1][j-1];
		}
	}
	for(j=1;j<lambda;j++){
		for(i=0;i<model->inputLen;i++){
			g[i][j] = e[0][j-1]*model->input[i] + g[i+1][j-1];
		}
	}
	//Ecuación (4.18)
	for(i=0;i<lambda;i++) h += g[0][i];
	//De modo similar, para fi, delta, mu, con alfa y beta.
	for(i=0;i<traject->outputLen;i++) fi[i][0] = traject->output[i];
	for(i=0;i<traject->inputLen;i++) delta [i][0] = traject->input[i];
	for(j=1;j<lambda;j++) fi[traject->outputLen][j-1] = 0;
	for(j=1;j<lambda;j++) delta[traject->inputLen][j-1] = 0;
	
	for(j=1;j<lambda;j++){
		for(i=0;i<traject->outputLen;i++){
			fi[i][j] = fi[0][j-1]*traject->output[i] + fi[i+1][j-1];
		}
	}
	for(j=1;j<lambda;j++){
		for(i=0;i<traject->inputLen;i++){
			delta[i][j] = fi[0][j-1]*traject->input[i] + delta[i+1][j-1];
		}
	}
	for(i=0;i<lambda;i++) mi += delta[0][i];
	
	//Ecuaciones B.19
	for(i=0; i<model->outputLen; i++){
		ni[i] = 0;
		for(j=0;j<lambda;j++) ni[i] += e[i][j];
	}
	for(i=0; i<model->inputLen; i++){
		gamma[i] = 0;
		for(j=0;j<lambda;j++) gamma[i] += g[i][j];
	}	
	
}

/*
*
*			MÉTODOS SET
*
*/


/*
*			Ajustar el set-point
*/
void APControl::setSetPoint(float sp){
  ySP[0] = sp;
}

/*
*			Ajustar la salida del proceso
*/
void APControl::setProcessOutput(float out){
  y[0] = out;
}

/*
*			Ajustar la acción de control manualmente
*/
void APControl::setControlAction(float action){
  u[0] = action;
}

void APControl::setLambda(byte l){
	lambda = l;
}

void APControl::setMaxControl(float m){
		maxControl = m;
}
	
void APControl::setMinControl(float m){
		minControl = m;
}
void APControl::setMaxIncControl(float m){
		maxIncControl = m;
}
void APControl::setDeltaB(float db){
		deltaB = db;
}
void APControl::setAdaptFact(float adapt){
		adaptFact = adapt;
}
void APControl::setModel(transferFunction *mod){
	model = mod;
}
void APControl::setTraject(transferFunction *traj){
	traject = traj;
}

/*
*
*			MÉTODOS GET
*
*/

float APControl::getYd(){
  return yDNext;
}

float APControl::getProcessOutput(){
  return y[0];
}

float APControl::getControlAction(){
  return u[0];
}

float APControl::getPostErr(){
  return err;
}

float APControl::getSetPoint(){
  return ySP[0];
}

int APControl::getLambda(){
	return lambda;
}

#ifdef DEBUG_MODE

float APControl::getAm1(){
	return model->output[0];
}

float APControl::getAm2(){
	return model->output[1];
}
float APControl::getAm3(){
	return model->output[2];
}
float APControl::getBm1(){
	return model->input[0];
}
float APControl::getBm2(){
	return model->input[1];
}
float APControl::getBm3(){
	return model->input[2];
}

#endif

