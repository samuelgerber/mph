//author: Samuel Gerber

#ifndef FILTRATION_H
#define FILTRATION_H


#include <map>
#include "Simplex.h"

//NOT IN USE

//Note vertexTime could be removed if 0-simplicies are added to the filtration.
//For computational speed this is not done here. Would require to adapt all the
//Rips, MutiscaleRips and PersistentHomology code if vertexTime is removed

typedef std::map<Simplex, double> IFiltration;
typedef IFiltration::iterator IFiltrationIterator;

class FiltrationEntry{
  public:
    Simplex simplex;
    double time;


    FiltrationEntry(const Simplex &s, double t) : simplex(s), time(t){
    }; 

    bool operator == (const FiltrationEntry& other) const{
      return time == other.time && simplex == other.simplex; 
    };

    bool operator < (const FiltrationEntry& other) const{
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

    bool operator > (const FiltrationEntry& other) const{
      return other < *this;
    };
};


typedef std::vector< FiltrationEntry  > Filtration;
typedef Filtration::iterator FiltrationIterator;





#endif 

