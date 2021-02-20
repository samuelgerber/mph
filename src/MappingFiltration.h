//author: Samuel Gerber

#ifndef MAPPINGFILTRATION_H
#define MAPPINGFILTRATION_H


#include <map>
#include <vector>
#include "MappingSimplex.h"


//Note vertexTime could be removed if 0-simplicies are added to the filtration.
//For computational speed this is not done here. Would require to adapt all the
//Rips, MutiscaleRips and PersistentHomology code if vertexTime is removed

typedef std::map<MappingSimplex, double> IMappingFiltration;
typedef IMappingFiltration::iterator IMappingFiltrationIterator;

class MappingFiltrationEntry{
  public:
    MappingSimplex simplex;
    double time;

    MappingFiltrationEntry(){
    };

    MappingFiltrationEntry(const MappingSimplex &s, double t):simplex(s), time(t){
    }; 

    bool operator == (const MappingFiltrationEntry& other) const{
      return time== other.time && simplex == other.simplex; 
    };

    bool operator < (const MappingFiltrationEntry& other) const{
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

    bool operator > (const MappingFiltrationEntry& other) const{
      return other < *this;
    };
};


typedef std::vector< MappingFiltrationEntry  > MappingFiltration;
typedef MappingFiltration::iterator MappingFiltrationIterator;





#endif 

