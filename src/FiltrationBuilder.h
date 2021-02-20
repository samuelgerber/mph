//author: Samuel Gerber

#ifndef FILTRATIONBUILDER_H
#define FILTRATIONBUILDER_H


#include <vector>
#include <limits>
#include <map>

#include "Filtration.h"



template<typename TPrecision>
class FiltrationBuilder{
  public:
    typedef typename std::map<int, TPrecision> NeighborMap;
    typedef typename NeighborMap::iterator NeighborMapIterator;
    typedef typename std::vector<NeighborMap> Neighbors; 
    

  private:
    //Simplex and adding time
    Filtration filtration;

  public:


    FiltrationBuilder(){ 

    };


    ~FiltrationBuilder(){

    };



    void run( Neighbors &N, 
              int maxD, 
              TPrecision *vertexTimes = NULL, 
              TPrecision minTime=0, 
              TPrecision maxTime=std::numeric_limits<TPrecision>::max(),  
              int maxSize = std::numeric_limits<int>::max() ){

#ifdef VERBOSE
      clock_t t0 = clock();
#endif
      
      filtration.clear();

      //Build rips based on nn graph
      for(unsigned int i=0; i < N.size(); i++){
        Simplex simplex;
        simplex.vertices.insert(i);
        TPrecision t=0;
        if(vertexTimes != NULL){
          t = vertexTimes[i];
        }
        FiltrationEntry fe(simplex, t);

        NeighborMap &knn = N[i];
        NeighborMap lower;
        lowerNeighbors(i, knn, lower);
        addCofaces(fe, N, lower, maxD, minTime, maxTime);    
      }


#ifdef VERBOSE
      clock_t t1 = clock();
#endif

      std::sort(filtration.begin(), filtration.end() );

#ifdef VERBOSE
      clock_t t2 = clock();
#endif
      //MappingFiltrationIterator it = std::unique (filtration.begin(), filtration.end());
      //filtration.erase( it, filtration.end() );


#ifdef VERBOSE
      clock_t t3 = clock();
      std::cout << " rips inc. duplicates size: " << filtration.size() << std::endl;
      std::cout << " rips unique size: " << filtration.size() << std::endl;

      std::cout << " build time: " << ((double) t1-t0) / CLOCKS_PER_SEC << std::endl;
      std::cout << " sort time: " << ((double) t2-t1) / CLOCKS_PER_SEC << std::endl;
      std::cout << " unique time: " << ((double) t3-t2) / CLOCKS_PER_SEC << std::endl;
#endif
     /*
      if(filtration.size() > maxSize + size){
        MappingFiltrationIterator fIt  = filtration.begin();
        std::advance(fIt,maxSize + size);
        filtration.erase( fIt, filtration.end() );
        truncated = true;
      }
      */


    };



    Filtration &getFiltration(){
      return filtration;
    };

    TPrecision getMaxTime(){
      return filtration.back().time;
    };




      
 private:


 

    void lowerNeighbors(int index, NeighborMap &nn, NeighborMap &lower){
      lower.clear();
      for(NeighborMapIterator it = nn.begin(); it != nn.end(); ++it){
        if(it->first < index){
          lower.insert(*it);
        }
      } 
    };




    void addCofaces(FiltrationEntry &f, Neighbors &N, NeighborMap &nn, unsigned int maxD,
        double minTime, double maxTime){

      if(f.time > maxTime){
        return;
      }
      
      filtration.push_back(f);
      
      if(f.simplex.vertices.size() > (maxD+1) ){
        return;
      }


      for(NeighborMapIterator it = nn.begin(); it != nn.end(); ++it){

        //Create new simplex a new entrance time
        Simplex s = f.simplex;
        double time = f.time;
        NeighborMap &knn = N[it->first];
        for(std::set<int>::iterator vit = f.simplex.vertices.begin(); vit !=
            f.simplex.vertices.end(); ++vit){
          double tmp = knn[*vit];
          if(tmp > time){
            time = tmp;
          }
        }
        if(time < minTime){
          time = minTime;
        }

        s.vertices.insert(it->first);
        FiltrationEntry fe(s, time);

        //Compute lower neighbors intersections
        NeighborMap lower;
        lowerNeighbors(it->first, knn, lower);
        NeighborMap res;
        std::set_intersection( lower.begin(), lower.end(), nn.begin(), nn.end(),
                               std::inserter( res, res.begin()), lower.value_comp() );


        addCofaces(fe, N, res, maxD, minTime, maxTime);
      }

    };


};

#endif 

