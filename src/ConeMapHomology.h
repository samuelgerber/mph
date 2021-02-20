//author: Samuel Gerber

#ifndef CONEMAPHOMOLOGY_H
#define CONEMAPHOMOLOGY_H


#include <map>
#include <set>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <stdlib.h>

#include "Filtration.h"




template <typename TPrecision>
class ConeMapHomology{

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;




    class Event{
      public:
        Event(const TPrecision &b, const TPrecision &d, int dim, int s) : birth(b),
        death(d), dimension(dim), scale(s){};

        TPrecision birth;
        TPrecision death;
        int dimension;
        int scale;
    };

    typedef std::vector< Event > Events;
    typedef typename Events::iterator EventsIterator;


    typedef std::map<Simplex, int > IDMap;
    typedef typename IDMap::iterator IDMapIterator;

    typedef std::map<int, FiltrationEntry> FiltrationEntryMap;
    typedef typename FiltrationEntryMap::iterator FiltrationEntryMapIterator;

    class Column{
      public:
        Column(){
          time = 0;
        }
        std::list<int> entries;
        TPrecision time;
        int dim;
    };



  private:

    std::vector< Column > columns;
    std::vector< int > low2id;
    IDMap ids;
    TPrecision maximalEntranceTime;

    int nDuplicates;
    Events events;



    void getColumn(const Simplex &s, std::list<int> &col){

      std::list<Simplex> faces = s.getFaces();
      for(std::list<Simplex>::iterator it = faces.begin(); it != faces.end();
          ++it){

        if(ids[*it] == 0){

#ifdef VERBOSE
          std::cout << "moo" << std::endl;
          for(std::set<int>::iterator vit = s.vertices.begin(); vit != s.vertices.end(); ++vit){
            std::cout << *vit << std::endl;
          }
          std::cout << "moo" << std::endl;
          for(std::set<int>::iterator vit = (*it).vertices.begin(); vit != (*it).vertices.end(); ++vit){
            std::cout << *vit << std::endl;
          }
          std::cout << "moo" << std::endl;

          exit(0); //TODO
#endif

        }
        else{
          col.push_back( getID(*it) );
        }
      }


      col.sort();
     // std::sort(col.begin(), col.end());
    };




    int getLow(std::list<int> &col){
      if( col.empty() ){
        return 0;
      }
      return col.back();
    };




    int getID(const Simplex &s){
      IDMapIterator it = ids.find(s);
      if(it == ids.end()){
        int id = ids.size() + 1;
        ids[s] = id;
        return id;
      }
      return it->second;
    };






  public:

    ConeMapHomology(){
      nDuplicates = 0;
      columns.resize(1);
      low2id.resize(1);
    };



    //The filtration should not contain any vertices but start at the edges
    //This methods expects and "inverted filtration", i.e. each simplex mapping
    //to a time, the sorting is done as a preprocessing step here
    // (Note: this conversion might not be necessary since the simplicies in the
    // set are in lexicographical order)
    void run( IFiltration &ifilt, unsigned int maxDim, int scale){

#ifdef VERBOSE
      std::cout << "IFilt size: " << ifilt.size() << std::endl;
#endif
      Filtration filt;
      for(IFiltrationIterator it = ifilt.begin(); it != ifilt.end();
          ++it){
        filt.push_back( FiltrationEntry(it->first, it->second) );
      }
      std::sort(filt.begin(), filt.end() );

      run(filt,  maxDim, scale);
    };



    //The filtration should not contain any vertices but start at the edges
    void run(Filtration &filt, unsigned int maxDim, int scale){

      int filtSize = filt.size();
      columns.resize( columns.size() + filtSize );
      low2id.resize( low2id.size() + filtSize );

      //Do matrix reduction
      for(FiltrationIterator it = filt.begin(); it != filt.end(); ++it){

        Simplex &current = it->simplex;

        int id = getID(current);
        Column &c1 = columns[id];
        //Already added column
        if( !c1.entries.empty() ){
           nDuplicates++;
           continue;
        }
        //c1.simplex = current;
        c1.time = it->time;
        c1.dim = current.vertices.size();
        maximalEntranceTime = c1.time;


        //Do reduction with column c1
        //clock_t t1 =  clock();
        getColumn(current, c1.entries);
        int low = getLow( c1.entries );
        bool collision = low != 0;
        while( collision ){

          int lowInd = low2id[low];
          Column &c2 = columns[ lowInd ];

          symmetric_diff_inplace(c1.entries, c2.entries);

          low = getLow(c1.entries);
          collision = ( !c2.entries.empty() && low != 0 );
        }
        //Check possible cases of deaths
        if( current.vertices.size() > 1 && low > 0 ){ //Death
          low2id[ low ] = id;
          Column &cLow = columns[low];
          if( cLow.time < it->time ){
            events.push_back( Event( cLow.time, it->time,
                                     current.vertices.size() -2 , scale) );
          }

        }

      };



#ifdef VERBOSE
      std::cout << "#Events: " << events.size() << std::endl;
      std::cout << "Filtration size: " << columns.size() << std::endl;
      std::cout << "Duplicates: " << nDuplicates << std::endl;
#endif
    };



    Events &getEvents(){
      return events;
    };



    MatrixXp getDiagram(){
      MatrixXp ph(4, events.size());
      int index = 0;
      for(typename Events::iterator it = events.begin(); it != events.end();
          ++it, ++index){
        Event &e = (*it);
        //ph(0, index) =  e.s.scale;
        ph(0, index) =   e.birth;
        ph(1, index) =   e.death;
        ph(2, index) =   e.scale;
        ph(3, index) =   e.dimension;
      }
      return ph;
    };

    TPrecision getMaximalTime(){
      return maximalEntranceTime;
    };


    int getFiltrationSize(){
       return columns.size() - nDuplicates;
    };

  private:

    void symmetric_diff_inplace(std::list<int> &c1, std::list<int> &c2){
      std::list<int>::iterator it1 = c1.begin();
      std::list<int>::iterator it2 = c2.begin();

      while (true)
      {
        if( it1 == c1.end() ){
          c1.insert( it1, it2, c2.end() );
          return;
        }
        if( it2 == c2.end() ){
          return;
        }

        if( *it1 < *it2 ) {
          ++it1;
        }
        else if( *it2 < *it1 ) {
          c1.insert( it1, *it2 );
          ++it2;
        }
        else {
          it1 = c1.erase( it1 );
          ++it2;
        }
      }
    };



};

#endif

