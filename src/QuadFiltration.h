//author: Samuel Gerber

#ifndef QUADFILTRATION_H
#define QUADFILTRATION_H


#include <map>
#include "QuadSimplex.h"

//NOT IN USE

//Note vertexTime could be removed if 0-simplicies are added to the filtration.
//For computational speed this is not done here. Would require to adapt all the
//Rips, MutiscaleRips and PersistentHomology code if vertexTime is removed

typedef std::map<QuadSimplex, double> IQuadFiltration;
typedef IQuadFiltration::iterator IQuadFiltrationIterator;

class QuadFiltrationEntry{
  public:
    QuadSimplex simplex;
    double time;


    QuadFiltrationEntry(const QuadSimplex &s, double t) : simplex(s), time(t){
    }; 

    bool operator == (const QuadFiltrationEntry& other) const{
      return time == other.time && simplex == other.simplex; 
    };

    bool operator < (const QuadFiltrationEntry& other) const{
      if(time < other.time){
        return true;
      }
      else if(time > other.time){
        return false;
      }
      else{
        return simplex < other.simplex;
      }
    };

    bool operator > (const QuadFiltrationEntry& other) const{
      return other < *this;
    };
};


typedef std::vector< QuadFiltrationEntry  > QuadFiltration;
typedef QuadFiltration::iterator QuadFiltrationIterator;





#endif 

