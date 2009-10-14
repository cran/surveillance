/**
 * File based on algoritmos.cpp and sv.cpp from the TerraView plugin.
 * C++ source originally created by Marcos Oliveira Prates on 06 April 2006
 *
 * R interface by Michael Höhle initiated on 12 Jan 2009
 */

#include "stcd-assuncaocorrea.h"
#include <cmath>
#include <iostream>

using namespace std;

// Calculate the number of events centered in xi,yi trough the critical
// radius given, from  the time tj until ti
int CalculaNCj(short **MSpace, const int EvtN, const int EvtJ)
{
  int i;
  int Soma=0;
  for (i=EvtJ;i<=EvtN;i++)
    Soma += MSpace[EvtJ][i];
  return(Soma);
}

// Calculate the number of events centered in xi,yi trough the critical
// radius given
int ContaEvt(short **MSpace, const int EvtN, const int EvtJ)
{
  int i;
  int Soma=0;
  for (i=0;i<=EvtN;i++)
    Soma += MSpace[EvtJ][i];
  return(Soma);
}

//////////////////////////////////////////////////////////////////////
// Unfortunately, this function has not been commented in the TerraView
// version
//////////////////////////////////////////////////////////////////////
int SistemadeVigilancia(SVEventLst &ev, const double RaioC, const double epslon,
                        std::valarray<double> &R) {
  size_t i, j, NCj, NumTotEvt, NumEvtCil;
  short **MSpace;
  double pontox, pontoy, DistEucl, Soma, UCj, fator;
 
  //order the event list
  ev.sort();
   
  SVEventLst::size_type n_event = ev.size();

  //create the spatio matrix
  MSpace = new short* [n_event];
  if( MSpace == NULL )
    return 1;
  for( i = 0; i < n_event; i++ ) {
    MSpace[i] = new short[n_event];
    if( MSpace[i] == NULL ) {
      delete []MSpace;
      return 1;
    }
  }

  //create the output vector
  R.resize(n_event);
  if( R.size() != n_event ) {
    for( i = 0; i < n_event; i++ ) {
      delete []MSpace[i];
    }
    delete []MSpace;
    return 1;
  }

  //populate the spatio matrix with 1 if is close in spatio and 0 if not  
  i = 0;
  for( SVEventLst::iterator it = ev.begin(); it != ev.end(); ++it, i++ ) {
    j = 0;
    for( SVEventLst::iterator jt = ev.begin(); jt != ev.end(); ++jt, j++ ) {
      pontox = (*it).x-(*jt).x;
      pontoy = (*it).y-(*jt).y;
      DistEucl = sqrt((pontox*pontox)+(pontoy*pontoy));
      if((DistEucl < RaioC))
	MSpace[i][j]=1;
      else
	MSpace[i][j]=0;
    }
  }

  //do the calculs to find the output value of each event
  for( i = 0; i < n_event; i++ ) {
    Soma = 0.0;
    for( j = 0; j <= i; j++ ) {
      NCj = CalculaNCj(MSpace,i,j);
      NumTotEvt = ContaEvt(MSpace,i,j);
      NumEvtCil = i-j+1;
      UCj = ((double)NumEvtCil*(double)NumTotEvt)/(double)(i+1);
      fator = 1.0+epslon;
      Soma += (pow(fator,(double)NCj) * exp((-epslon)*UCj));
    }
    R[i] = Soma;
  }

  //clean memory
  for( i = 0; i < n_event; i++ ) {
    delete [] MSpace[i];
  }
  delete [] MSpace;
  return 0;
}

int CalculaLambda(SVEventLst &ev, const double RaioC, const double epslon, std::valarray<double> &R, unsigned int &numObs)
{
  size_t i, j, NCj, NumTotEvt, NumEvtCil;
  short **MSpace;
  double pontox, pontoy, DistEucl, UCj, fator, lambda, lambdaMax;
  
  ev.sort();
  
  SVEventLst::size_type n_event = ev.size();
  
  //create the spatio matrix
  MSpace = new short* [n_event];
  if( MSpace == NULL )
    return 1;
  for( i = 0; i < n_event; i++ ) {
    MSpace[i] = new short[n_event];
    if( MSpace[i] == NULL ) {
      delete []MSpace;
      return 1;
    }
  }
  
  //create the output vector
  R.resize(n_event);
  if( R.size() != n_event ) {
    for( i = 0; i < n_event; i++ ) {
      delete []MSpace[i];
    }
    delete []MSpace;
    return 1;
  }
  
  //populate the spatio matrix with 1 if is close in spatio and 0 if not  
  i = 0;
  for( SVEventLst::iterator it = ev.begin(); it != ev.end(); ++it, i++ ) {
    j = 0;
    for( SVEventLst::iterator jt = ev.begin(); jt != ev.end(); ++jt, j++ ) {
      pontox = (*it).x-(*jt).x;
      pontoy = (*it).y-(*jt).y;
      DistEucl = sqrt((pontox*pontox)+(pontoy*pontoy));
      if((DistEucl < RaioC))
	MSpace[i][j]=1;
      else
	MSpace[i][j]=0;
    }
  }
  
  //do the calculus to find the output value of each event  
  i = numObs;
  lambdaMax = 0;
  for( j = 0; j <= i; j++ ) {
    NCj = CalculaNCj(MSpace,i,j);
    NumTotEvt = ContaEvt(MSpace,i,j);
    NumEvtCil = i-j+1;
    UCj = ((double)NumEvtCil*(double)NumTotEvt)/(double)(i+1);
    fator = 1.0+epslon;
    lambda = (pow(fator,(double)NCj) * exp((-epslon)*UCj));
    if (lambda > lambdaMax){
      lambdaMax = lambda;
      numObs = j;
    }
  }
  
  
  
  //clean memory
  for( i = 0; i < n_event; i++ ) {
    delete [] MSpace[i];
  }
  delete [] MSpace;
  return 0;
}

//////////////////////////////////////////////////////////////////////
// Shiryaev-Roberts space time detection as explained in the paper
// by Correa and Assuncao (2009).
//
// Params:
//  x - array with x location of events
//  y - array with y location of events
//  t - array with time point of the events (on some arbitrary time scale)
//  n - number of elements in x, y and t (the same for the three vectors)
//  radius  - cluster of the radius
//  epsilon - relative ratio of the intensity functions to detect for
//  threshold -- upper threshold when to sound the alarm
//  Rarray -- array of length n, this will contain the statistics calced
//            by the function
// idxFirstAlarm -- index in the x,y,t vector resulting in the alarm
// idxClusterCenter -- index in the x,y,t vector containing the cluster
//                      center
//////////////////////////////////////////////////////////////////////


extern "C" {

  void SRspacetime(double *x, double *y, double *t, long *n, double *radius,
		   double *epsilon, double *threshold, 
		   double *Rarray, int *idxFirstAlarm, int *idxClusterCenter) {
  
  //Create SVEventLst
  SVEvent e;
  SVEventLst eList;
  unsigned int i;
  long j;
  //Fill coordinates of event list
  for(j=0;j<*n;j++){
    e.x = x[j];
    e.y = y[j];
    e.t = t[j];
    eList.push_back(e);
  }
	
  //Array of test statistic values
  std::valarray<double> R;
  
  //Call SistemadeVigilancia
  SistemadeVigilancia(eList,*radius,*epsilon,R);

  //Debug purposes
  //cout << "Size of R = " << R.size() << endl;
  //Move values of test statistic for return 
  for(i=0;i<R.size();i++){
    //cout << "R[" << i << "] = " << R[i] << endl;
    Rarray[i] = R[i];
  }
  
  //Find index, where R is above the limit for the first time
  //Boolean "controle" indicates whether there is an alarm or not.
  bool controle = false;
  for(i=0;i<R.size();i++){
    if(R[i]>*threshold){
      controle = true;
      break;
    }
  }

  //Advancing the iterator it to the point
  //where the alarm is generated. 
  if (controle) {
    unsigned int cont = 0;
    SVEventLst::iterator it = eList.begin();
    while((cont < i) &&  (it != eList.end())){
      ++it;
      ++cont;
    }
    *idxFirstAlarm = cont;

    //Determine the cluster center of the alarm
    unsigned int num = cont;	
    CalculaLambda(eList,*radius,*epsilon,R,num);

    //Index of the cluster center
    *idxClusterCenter = num;
  } else {
    //If no alarms, then return -1 for both alarm idx and cluster center idx
    *idxFirstAlarm = -2;
    *idxClusterCenter = -2;
  }
  //Clean up (nothing to clean) and done
}

}
