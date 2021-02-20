//author: Samuel Gerber

#ifndef PERSISTENTHOMOLOGY_H
#define PERSISTENTHOMOLOGY_H


#include <map>
#include <set>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <stdlib.h>

#include "Simplex.h"
#include "Filtration.h"


#include <Eigen/Dense>

template <typename TPrecision>
class PersistentHomology{
 
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  
  
    class Event{
      public:
        Event(const Simplex &sIn, const TPrecision &b, const TPrecision &d) : s(sIn), birth(b),
        death(d){};

      Simplex s;
      TPrecision birth;
      TPrecision death;
    };

    typedef std::list< Event > Events;
    typedef typename Events::iterator EventsIterator;





  private:

      Events events;
      std::set< FiltrationEntry > elves;

      void getColumn(const Simplex &s, std::map<Simplex, int>
          &ids, std::list<int> &col){

        std::list<Simplex> faces = s.getFaces();
        for(std::list<Simplex>::iterator it = faces.begin(); it != faces.end();
            ++it){

          std::map<Simplex, int>::iterator fit = ids.find(*it);
          if( fit->second== 0){
#ifdef VERBOSE
            std::cout << "moo" << std::endl;
#endif
          }
          col.push_back( fit->second );
        }
        col.sort();
     };


      int getLow(std::list<int> &col){
        if(col.size() == 0){
          return 0;
        }
        return col.back();
      };



  public:
  
    PersistentHomology(){ 
    };
 


    //The filtration should not contain any vertices but start at the edges
    //This methods expects and "inverted filtration", i.e. each simplex mapping
    //to a time the sorting is done as a preprocessing step here
    // (Note: this conversion might not be necessary since the simplicies in the
    // set are in lexicographical order)
    void run(IFiltration &ifilt){
#ifdef VERBOS
      std::cout << "IFilt size: " << ifilt.size() << std::endl;
#endif

      clock_t start = clock();
      clock_t end = clock();
      unsigned long millis = (end - start) * 1000 / CLOCKS_PER_SEC;

      start = clock();
      Filtration filt;
      for(IFiltrationIterator it = ifilt.begin(); it != ifilt.end();
          ++it){
        filt.push_back( FiltrationEntry(it->first, it->second) );
      }
      std::sort( filt.begin(), filt.end() );
      //filt.sort();
      
      run( filt, filt.size() );

    };
    



    //The filtration should not contain any vertices but start at the edges
    //If maxD  > 0, "undead" events up to homology maxD are recored as well 
    void run( Filtration &filt, int maxD = 0){

      //Simplex to ID mapping
      std::map<Simplex, int> ids;

      //one indexing, zero reserved for empty column
      //Low mapping to Simplex ID
      std::vector< int > low2id( filt.size()+1, 0);

      //Index adding time
      std::vector<TPrecision> time( filt.size()+1, 0);

      //columns, id to a set of ids
      std::vector< std::list<int> > columns( filt.size()+1 );

      //Do matrix reduction
      events.clear();
      elves.clear();
      int nZeros = 0;
      TPrecision currentTime = 0;
      std::cout << filt.size() << std::endl;
      for(FiltrationIterator it = filt.begin(); it != filt.end(); ++it){

        if( it->time < currentTime){
#ifdef VERBOSE
          std::cout << "harumpf" << std::endl;
#endif
        }
        currentTime = it->time;

        const Simplex &current = it->simplex;

        std::list<int> cIds;
        getColumn(current, ids, cIds);
        int low = getLow(cIds);
        bool collision =  low != 0;
        while( collision ){
          std::list<int> &col = columns[ low2id[low] ];
          std::list<int> diff;
          std::set_symmetric_difference(cIds.begin(), cIds.end(), col.begin(),
              col.end(), std::inserter(diff, diff.end()) );
          cIds = diff;
          low = getLow(cIds);
          collision = (col.size() != 0 && low != 0) ;
        }


        // one indexing since 0 is reserved for empty column
        int index = ids.size() + 1;
        ids[ current ] = index;
        //id2low[ index ] = low;
        time[ index ] = it->time;
        columns[index] = cIds;
        if( low != 0 ){
          low2id[ low ] = index;
          if( it->time - time[low] > 0 ){
            events.push_back( Event( current, time[low], it->time ) );
          }
          if( maxD > 0){
            elves.erase( filt[ low-1] );
          }
        }
        else{

          if( it->simplex.vertices.size() > 1 && it->simplex.vertices.size() < maxD+2){
            elves.insert( *it ) ;
          }
          nZeros++;
        }

      };


#ifdef VERBOSE
      std::cout << "#Zeros: " << nZeros << std::endl;
      std::cout << "#Events: " << events.size() << std::endl;
#endif
    
    };



    Events &getEvents(){
      return events; 
    };



    MatrixXp getDiagram(){
      
      MatrixXp ph(3, events.size() + elves.size());
      
      int index = 0;
      for(typename Events::iterator it = events.begin(); it != events.end();
          ++it, ++index){
        Event &e = (*it);
        //ph(0, index) =  e.s.scale;
        ph(0, index) =  e.birth;
        ph(1, index) =  e.death;
        ph(2, index) =  e.s.vertices.size() - 2; 
      }
      for(typename std::set< FiltrationEntry >::iterator it = elves.begin(); it !=elves.end(); ++it, ++index){
        ph(0, index) =  it->time;
        ph(1, index) =  -1;
        ph(2, index) =  it->simplex.vertices.size() - 1; 
      }

      return ph;
    };



};

#endif 

