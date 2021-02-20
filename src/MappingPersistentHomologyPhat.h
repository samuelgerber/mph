//author: Samuel Gerber

#ifndef MAPPINGPERSISTENTHOMOLOGYPHAT_H
#define MAPPINGPERSISTENTHOMOLOGYPHAT_H


#include <map>
#include <set>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <stdlib.h>

#include "MappingFiltration.h"

//PHAT stuff

// wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm
#include <phat/compute_persistence_pairs.h>

// main data structure (choice affects performance)
#include <phat/representations/vector_vector.h>

// algorithm (choice affects performance)
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/twist_reduction.h>



template <typename TPrecision>
class MappingPersistentHomologyPhat{

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;


    typedef std::map<MappingSimplex, int > IDMap;
    typedef typename IDMap::iterator IDMapIterator;

    typedef std::map<int, MappingFiltrationEntry> FiltrationEntryMap;
    typedef typename FiltrationEntryMap::iterator FiltrationEntryMapIterator;



  private:

    MatrixXp dgm;
    IDMap ids;
    TPrecision maximalEntranceTime;
    FiltrationEntryMap active;


    void getColumn(const MappingSimplex &s, phat::column &col){

      std::list<MappingSimplex> faces = s.getFaces();
      for(std::list<MappingSimplex>::iterator it = faces.begin(); it != faces.end();
          ++it){
        col.push_back( getID(*it) );
      }

      std::sort(col.begin(), col.end());
    };



    int getID(const MappingSimplex &s){
      IDMapIterator it = ids.find(s);
      if(it == ids.end()){
        int id = ids.size();
        ids[s] = id;
        return id;
      }
      return it->second;
    };





  public:

    MappingPersistentHomologyPhat(){
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
    void run( MappingFiltration &filt, TPrecision mapTime, unsigned int maxDim,
              bool storeActie= false){

      //clear any possible leftovers
      ids.clear();

      int filtSize = filt.size();
      std::vector<double> times( filtSize );
      std::vector<int> scales( filtSize );
      std::vector<bool> mapped( filtSize );

      phat::boundary_matrix< phat::vector_vector > boundary_matrix;
      boundary_matrix.set_num_cols( filtSize );
      std::vector< phat::index > temp_col;
      for(MappingFiltrationIterator it = filt.begin(); it != filt.end(); ++it){
        MappingSimplex &current = it->simplex;
        int id = getID(current);
        times[id] = it->time;
        scales[id] = current.getScale();
        mapped[id]  = current.isMapped();
        int dim = current.vertices.size();
        boundary_matrix.set_dim(id, dim-1);
        if(dim == 1){
          boundary_matrix.set_col( id, temp_col);
        }
        else{
          getColumn(current, temp_col);
          boundary_matrix.set_col( id, temp_col);
          temp_col.clear();
        }
      }
      // define the object to hold the resulting persistence pairs
      phat::persistence_pairs pairs(times);

      // choose an algorithm (choice affects performance) and compute the persistence pair
      // (modifies boundary_matrix)
      phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );

      //Store persistence diagram in matrix
      int index = 0;
      dgm = MatrixXp(4, pairs.get_num_pairs() );
      for(int i=0; i<pairs.get_num_pairs(); i++){
        int birth = pairs.get_pair( i ).first;
        int death = pairs.get_pair( i ).second;

        if( !mapped[death] ){
          dgm(0, index) = times[birth];
          dgm(1, index) = times[death];
          dgm(3, index) = scales[death];
          dgm(2, index) = boundary_matrix.get_dim( birth );
          ++index;
        }
      }
      dgm = dgm.leftCols(index);

      maximalEntranceTime = times.back();

#ifdef VERBOSE
      std::cout << "#Active: " << active.size() << std::endl;
#endif
    };


    MatrixXp getDiagram(){
      return dgm;
    };

    TPrecision getMaximalTime(){
      return maximalEntranceTime;
    };

    FiltrationEntryMap &getActive(){
      return active;
    };





};

#endif

