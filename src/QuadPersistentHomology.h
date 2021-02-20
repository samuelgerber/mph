//author: Samuel Gerber

#ifndef QUADPERSISTENTHOMOLOGY_H
#define QUADPERSISTENTHOMOLOGY_H


#include <map>
#include <set>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <stdlib.h>
#include <math.h>

#include "QuadFiltration.h"


#include <Eigen/Dense>

template <typename TPrecision>
class QuadPersistentHomology{
 
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  
  
    class Event{
      public:
        Event(const QuadSimplex &bs, const QuadSimplex &ds, const TPrecision &b, const TPrecision &d) : 
           sb(bs), sd(ds), birth(b),  death(d){};

      QuadSimplex sb;
      QuadSimplex sd;
      TPrecision birth;
      TPrecision death;
    };

    typedef std::list< Event > Events;
    typedef typename Events::iterator EventsIterator;





  private:

      Events events;

      void getColumn(const QuadSimplex &s, std::map<QuadSimplex, int>
          &ids, std::list<int> &col){

        std::list<QuadSimplex> faces = s.getFaces();
        for(std::list<QuadSimplex>::iterator it = faces.begin(); it != faces.end();
            ++it){

          std::map<QuadSimplex, int>::iterator fit = ids.find(*it);
          //std::cout << "Face: " << it->vertices.size() << std::endl;
          if( fit != ids.end() ){
            if( fit->second== 0){
#ifdef VERBOSE
              std::cout << "moo" << std::endl;
#endif
            }
            //std::cout << s.vertices.size() << " : " << fit->second << std::endl;
            col.push_back( fit->second );
            
          }
        }
        //std::cout << " - " << std::endl;
        col.sort();
     };


      int getLow(std::list<int> &col){
        if(col.size() == 0){
          return 0;
        }
        return col.back();
      };



  public:
  
    QuadPersistentHomology(){ 
    };
 


    //The filtration should not contain any vertices but start at the edges
    //This methods expects and "inverted filtration", i.e. each simplex mapping
    //to a time the sorting is done as a preprocessing step here
    // (Note: this conversion might not be necessary since the simplicies in the
    // set are in lexicographical order)
    void run(IQuadFiltration &ifilt){
#ifdef VERBOS
      std::cout << "IFilt size: " << ifilt.size() << std::endl;
#endif

      clock_t start = clock();
      clock_t end = clock();
      unsigned long millis = (end - start) * 1000 / CLOCKS_PER_SEC;

      start = clock();
      QuadFiltration filt;
      for(IQuadFiltrationIterator it = ifilt.begin(); it != ifilt.end();
          ++it){
        filt.push_back( QuadFiltrationEntry(it->first, it->second) );
      }
      std::sort( filt.begin(), filt.end() );
      //filt.sort();
      
      run( filt, filt.size() );

    };
    



    //The filtration should not contain any vertices but start at the edges
    void run( QuadFiltration &filt ){

      //QuadSimplex to ID mapping
      std::map<QuadSimplex, int> ids;

      //one indexing, zero reserved for empty column
      //Low mapping to QuadSimplex ID
      std::vector< int > low2id( filt.size()+1, 0);

      //Index adding time
      std::vector<TPrecision> time( filt.size()+1, 0);

      //columns, id to a set of ids
      std::vector< std::list<int> > columns( filt.size()+1 );

      //Do matrix reduction
      int nZeros = 0;
      TPrecision currentTime = 0;
      std::cout << filt.size() << std::endl;
      for(QuadFiltrationIterator it = filt.begin(); it != filt.end(); ++it){

        if( it->time < currentTime){
#ifdef VERBOSE
          std::cout << "harumpf" << std::endl;
#endif
        }
        currentTime = it->time;

        const QuadSimplex &current = it->simplex;

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
            events.push_back( Event( filt[low-1].simplex, current, time[low], it->time ) );
          }
        }
        else{
          nZeros++;
        }

      };


#ifdef VERBOSE
      std::cout << "#zeros: " << nZeros << std::endl;
      std::cout << "#Events: " << events.size() << std::endl;
#endif
    
    };



    Events &getEvents(){
      return events; 
    };



    MatrixXp getDiagram(){
      
      MatrixXp ph(5, events.size());
      
      int index = 0;
      for(typename Events::iterator it = events.begin(); it != events.end();
          ++it, ++index){
        Event &e = (*it);
        //ph(0, index) =  e.s.scale;
        ph(0, index) =  e.birth;
        ph(1, index) =  e.death;
        ph(2, index) =  log2( e.sd.vertices.size() ) - 1;
        ph(3, index) =  *e.sb.vertices.begin(); 
        ph(4, index) =  *e.sd.vertices.begin(); 
      }

      return ph;
    };



};

#endif 

