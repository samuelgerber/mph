//author: Samuel Gerber

#ifndef MAPPINGPERSISTENTHOMOLOGY_H
#define MAPPINGPERSISTENTHOMOLOGY_H


#include <map>
#include <set>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <stdlib.h>

#include "MappingFiltration.h"




template <typename TPrecision>
class MappingPersistentHomology{

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;

    


    class Event{
      public:
        Event(const TPrecision &b, const TPrecision &d, int dim) : birth(b),
        death(d), dimension(dim){};

        TPrecision birth;
        TPrecision death;
        int dimension;
    };

    typedef std::vector< Event > Events;
    typedef typename Events::iterator EventsIterator;


    typedef std::map<MappingSimplex, int > IDMap;
    typedef typename IDMap::iterator IDMapIterator;

    typedef std::map<int, MappingFiltrationEntry> FiltrationEntryMap;
    typedef typename FiltrationEntryMap::iterator FiltrationEntryMapIterator;
    
    class Column{
      public:
        Column(){
          time = 0;
        }
        std::list<int> entries;
        TPrecision time;
        int dim;
        //MappingSimplex simplex;
    };



  private:

    IDMap ids;
    IDMap mids;
    TPrecision maximalEntranceTime;

    Events events;

    FiltrationEntryMap active;


    void getColumn(const MappingSimplex &s, std::list<int> &col, int filtSize){

      std::list<MappingSimplex> faces = s.getFaces();
      for(std::list<MappingSimplex>::iterator it = faces.begin(); it != faces.end();
          ++it){
        if(ids[*it] == 0){
#ifdef VERBOSE
          std::cout << "moo" << std::endl;
#endif
        }
        col.push_back( getID(*it) );
      }

      //Add mapped simplex id at top of column with negative id
      /*
      MappingSimplex mapped;
      mapped.vertices = s.getMapsTo();
      if( mapped.vertices.size() == s.vertices.size() ){
        col.push_back( getMappedID(mapped) - filtSize - 1 );
      }
      */
      
      col.sort();
     // std::sort(col.begin(), col.end());
    };




    int getLow(std::list<int> &col){
      if( col.empty() ){
        return 0;
      }
      return col.back();
    };




    int getID(const MappingSimplex &s){
      IDMapIterator it = ids.find(s);
      if(it == ids.end()){
        int id = ids.size() + 1;
        ids[s] = id;
        return id;
      }
      return it->second;
    };



/*
    int getMappedID(const MappingSimplex &s){
      IDMapIterator it = mids.find(s);
      if(it == mids.end()){
        int id = mids.size() + 1;
        mids[s] = id;
        return id;
      }
      return it->second;
    };

*/



  public:

    MappingPersistentHomology(){ 
    };



    //The filtration should not contain any vertices but start at the edges
    //This methods expects and "inverted filtration", i.e. each simplex mapping
    //to a time, the sorting is done as a preprocessing step here
    // (Note: this conversion might not be necessary since the simplicies in the
    // set are in lexicographical order)
    void run(IMappingFiltration &ifilt, TPrecision mapTime, unsigned int maxDim){

#ifdef VERBOSE
      std::cout << "IFilt size: " << ifilt.size() << std::endl;
#endif
      MappingFiltration filt;
      for(IMappingFiltrationIterator it = ifilt.begin(); it != ifilt.end();
          ++it){
        filt.push_back( MappingFiltrationEntry(it->first, it->second) );
      }
      std::sort(filt.begin(), filt.end() );
      //filt.sort();

      run(filt, mapTime, maxDim);
    };



    //The filtration should not contain any vertices but start at the edges
    void run(MappingFiltration &filt, TPrecision mapTime, unsigned int maxDim){

      //clear any possible leftovers
      ids.clear();
      mids.clear();
      events.clear();
      active.clear();

      int filtSize = filt.size();
      //columns, id to a set of ids
      std::vector< Column > columns( filtSize + 1 );
      std::vector< int > low2id( filtSize + 1 );
      //std::vector< int > mLow2id( filtSize + 1 );
      
      

      //Do matrix reduction
      TPrecision currentTime = 0;
      //FortranLinalg::DenseVector<int> delta(filtSize);
      for(MappingFiltrationIterator it = filt.begin(); it != filt.end(); ++it){

        MappingSimplex &current = it->simplex;
        int id = getID(current);
        Column &c1 = columns[id];
        //c1.simplex = current;
        c1.time = it->time;
        c1.dim = current.vertices.size();
        currentTime = c1.time;


        //Do reduction with column c1
        //clock_t t1 =  clock();
        getColumn(current, c1.entries, filtSize);
        int low = getLow( c1.entries );
        bool collision = low != 0;
        while( collision ){
         
          int lowInd = 0;
          //if(low < 0){
          //  lowInd = mLow2id[-low];
          //}
          //else{
          lowInd = low2id[low];
          //} 
          Column &c2 = columns[ lowInd ];
          

          /*
          std::list<int> diff;
          std::set_symmetric_difference(c1.entries.begin(), c1.entries.end(),
             c2.entries.begin(), c2.entries.end(), std::inserter(diff, diff.end()) );
          c1.entries = diff;
          */
          symmetric_diff_inplace(c1.entries, c2.entries);

          low = getLow(c1.entries);
          collision = ( !c2.entries.empty() && low != 0 );
        }
        //clock_t t2 =  clock();
        //delta(index) = (t2-t1);



        //Check possible cases of deaths
        if(current.vertices.size() > 1){
          if( low > 0 ){ //Death
            
            low2id[ low ] = id;
            Column &cLow = columns[low];
            if( cLow.time < it->time  && !current.isMapped() ){//!maps[id] ){
              events.push_back( Event( cLow.time, it->time,
                    current.vertices.size()-2 ) );
            }
            //else{
             // active.erase( low );
            //}
          
          }
          else{//Birth
            if(current.vertices.size() < (maxDim+2) ){
              //Potentially dies in mapping
              //current.setActive( true );
              //if( !current.isMapped() ){
              //active[id] = *it;
              //}
              //activeCount[id] = count;
              //active.insert(id);
              //std::cout << current.vertices.size() << std::endl;
              //std::cout << "#Active: " << zeroes.size() << std::endl;
            }
            //Column &c1 = columns[id];
            //if(c1.time < mapTime && current.vertices.size() < (maxDim+2) ){
              //Potentially dies in mapping
              //zeroes.insert(id);
            //}
          }
          //else{
            //it will be either transfered or die within this scale
          //  mLow2id[-low] = id;
          //}
        }
      
      };


      maximalEntranceTime = currentTime;
      
      //Detect deaths in mapping - Can't happen?
      /*
      int nDM = 0;
      for(std::set<int>::iterator it = zeroes.begin(); it != zeroes.end();
          ++it){
        int id = *it;
        Column &c1 = columns[id];
        if(c1.time < mapTime ){
          nDM++;
          events.push_back( Event( c1.time, mapTime, c1.dim - 1  ) );
        }
      }
      */


      //FortranLinalg::LinalgIO<int>::writeVector("times.data", delta);
      //std::cout << "#zeros: " << zeroes.size() << std::endl;
#ifdef VERBOSE
      std::cout << "#Events: " << events.size() << std::endl;
      std::cout << "#Active: " << active.size() << std::endl;
#endif
    };



    Events &getEvents(){
      return events; 
    };



    MatrixXp getDiagram(){
      MatrixXp ph(3, events.size());
      int index = 0;
      for(typename Events::iterator it = events.begin(); it != events.end();
          ++it, ++index){
        Event &e = (*it);
        //ph(0, index) =  e.s.scale;
        ph(0, index) =   e.birth;
        ph(1, index) =   e.death;
        ph(2, index) =   e.dimension;
      }
      return ph;
    };

    TPrecision getMaximalTime(){
      return maximalEntranceTime;
    };

    FiltrationEntryMap &getActive(){
      return active;
    };



  private:

    void symmetric_diff_inplace(std::list<int> &c1, std::list<int> &c2){
      std::list<int>::iterator it1 = c1.begin();
      std::list<int>::iterator it2 = c2.begin();

      while (true)
      {
        if( it1 == c1.end() ){
          c1.insert(it1, it2, c2.end() );
          return;
        }
        if( it2 == c2.end() ){
          return;
        }

        if( *it1 < *it2 ) { 
          ++it1; 
        }
        else if( *it2 < *it1 ) { 
          c1.insert(it1, *it2);
          ++it2;
        }
        else { 
          it1 = c1.erase(it1); 
          ++it2; 
        }
      }
    };



};

#endif 

