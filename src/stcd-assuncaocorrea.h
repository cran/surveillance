/**
 * File based on algoritmos.cpp and sv.cpp from the TerraView plugin.
 * C++ source originally created by Marcos Oliveira Prates from the
 * Department of Statistics, UFMG, Brazil on 06 April 2006
 * 
 * R interface by Michael Höhle initiated on 12 Jan 2009
 * Note: Some function names and documentation are in Portugese
 */


#ifndef SRSPACETIME_H
#define SRSPACETIME_H

#include <list>
#include <valarray>

struct SVEvent {
 double x, y, t;
 friend bool operator<(const SVEvent &a, const SVEvent &b) {
  return (a.t < b.t);
 }
};

//STL is used (check its use)
typedef std::list<SVEvent> SVEventLst;

//Functions provided in sr-spacetime.cc
int CalculaNCj(short **MSpace, const int EvtN, const int EvtJ);
int ContaEvt(short **MSpace, const int EvtN, const int EvtJ);
//int SistemadeVigilancia(SVEventLst &, const double RaioC, const double epslon, 
//			std::valarray<double> &R);
//New version with different estimation approach
int SistemadeVigilancia(SVEventLst &ev,
			const double RaioC, const double epslon,  
			const double areaA, double *areaAcapBk, 
			const int cusum, std::valarray<double> &R);
int CalculaLambda(SVEventLst &ev, const double RaioC, const double epslon, 
		  std::valarray<double> &R, unsigned int &numObs);



// Hoehle wrapper function to create SVEvent list
//void SRspacetime(double *x, double *y, double *t, int *n, double *radius, double *epsilon, double *Rarray);


#endif
