//author: Samuel Gerber

#ifndef CONEMAPFILTRATIONBUILDER2_H
#define CONEMAPFILTRATIONBUILDER2_H


#include <map>
#include <set>
#include <vector>
#include <utility>
#include <limits>

#include "Filtration.h"




//Computes a Rips filtration F2 given a metric space X2, a filtration F1 on
//X1 and a mapping g from X1 to X2. The filtration is computed by first
//inserting the simplicies from F1 that map to F2 according to their time in
//filtration F1 and then adding according to Rips on F2 any remaining simplicies
//S with times max( max( time( F1 ), time(S) ) )

template<typename TPrecision>
class ConeMapFiltrationBuilder2{


  public:
    typedef typename std::map<int, TPrecision> NeighborMap;
    typedef typename NeighborMap::iterator NeighborMapIterator;
    typedef typename std::map<int, NeighborMap> Neighbors;
    typedef typename Neighbors::iterator NeighborsIterator;



  private:
    IFiltration &cone;
    bool truncated;

  public:


    ConeMapFiltrationBuilder2(IFiltration &c) : cone(c) {
      truncated = false;
    };


    ~ConeMapFiltrationBuilder2(){
    };



    //Compute the filtration based on the nearest neighbor structure on X2, and
    //the minimal insertion time minTime = max(F1). Compute the Rips up to time
    //maxTime and with simplicial dimension up to maxD. Truncate the Rips at
    //maxSize. Annotate the newly added simplicies, i.e. not mapped
    //from F1, as coming from MappingFiltration at scale "scale". This information is
    //required to run a multicsale persistent homology.  Before this is run
    //setMappedMappingFiltration (setting F1) should be called. If no mapped filtration
    //is set this is simply doing Rips with using a minimal entrance time of
    //minTime, i.e. anything added before that will be changed to entrance time
    //minTime, but adds some additional overhead computations.
    void run( Neighbors &N, TPrecision minTime, TPrecision maxTime,
              unsigned int maxSize, int maxD, int scale, Filtration &filtration){

#ifdef VERBOSE_DETAIL
      clock_t t0 = clock();
#endif

      //Build rips based on nn graph and mapped filtration
      filtration.clear();
      filtration.reserve(N.size()*(maxD+2));
      for( NeighborsIterator it = N.begin(); it !=N.end(); ++it ){
        Simplex simplex;
        simplex.vertices.insert( it->first );
        FiltrationEntry fe(simplex, minTime);

        NeighborMap &knn = it->second;
        NeighborMap lower;
        lowerNeighbors(it->first, knn, lower );
        addCofaces(fe, N, lower, maxD, minTime, maxTime, filtration);
      }

#ifdef VERBOSE_DETAIL
      clock_t t1 = clock();
#endif

      std::sort( filtration.begin(), filtration.end() );

#ifdef VERBOSE_DETAIL
      clock_t t2 = clock();
      //clock_t t3 = clock();
      std::cout << " rips size: " << filtration.size() << std::endl;
      std::cout << " build time: "  << ((double) t1-t0) / CLOCKS_PER_SEC << std::endl;
      std::cout << " sort time: "   << ((double) t2-t1) / CLOCKS_PER_SEC << std::endl;
      //std::cout << " unique time: " << ((double) t3-t2) / CLOCKS_PER_SEC << std::endl;
#endif



    };


    bool isTruncated(){
      return truncated;
    };




  private:



    void lowerNeighbors(int index, NeighborMap &nn, NeighborMap &lower ){
      lower.clear();
      for(NeighborMapIterator it = nn.begin(); it != nn.end(); ++it){
        if(it->first < index){
          lower.insert( *it );
        }
      }
    };




    void addCofaces( FiltrationEntry &f, Neighbors &N, NeighborMap &nn,
                     unsigned int maxD, double minTime, double maxTime, Filtration &filtration){

      if(f.time > maxTime){
        return;
      }

      //if( cone.find(f.simplex) == cone.end() ){
        filtration.push_back(f);
      //}


      if(f.simplex.vertices.size() > (maxD+1) ){
        return;
      }


      for(NeighborMapIterator it = nn.begin(); it != nn.end(); ++it){

        //Create new simplex a new entrance time
        Simplex s = f.simplex;
        double time = f.time;
        NeighborMap &knn = N[ it->first ];
        for(std::set<int>::iterator vit = f.simplex.vertices.begin(); vit !=
            f.simplex.vertices.end(); ++vit){
          double tmp = knn[ *vit ];
          if(tmp > time){
            time = tmp;
          }
        }
        if(time < minTime){
          time = minTime;
        }

        s.vertices.insert( it->first );
        FiltrationEntry fe(s, time);

        //Compute lower neighbors intersections
        NeighborMap lower;
        lowerNeighbors( it->first, knn, lower);
        NeighborMap res;
        std::set_intersection( lower.begin(), lower.end(), nn.begin(), nn.end(),
                               std::inserter( res, res.begin()), lower.value_comp() );


        addCofaces(fe, N, res, maxD, minTime, maxTime, filtration );

      }

    };

};

#endif

