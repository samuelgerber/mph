//author: Samuel Gerber


#ifndef CONEMAPRIPSPHAT_H
#define CONEMAPRIPSPHAT_H


#include "GMRATree.h"
#include "GMRANeighborhood.h"
#include "RelativeGMRANeighborhood.h"
#include "GMRAIDNode.h"

#include "ConeMapFiltrationBuilder2.h"

#include <map>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <vector>

#ifdef SAVETOFILE
#include <Eigen/Dense>
#include "EigenLinalgIO.h"
#endif

//PHAT stuff

// wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm
#include <phat/compute_persistence_pairs.h>

// main data structure (choice affects performance)
#include <phat/representations/vector_vector.h>
#include <phat/representations/vector_heap.h>
#include <phat/representations/vector_set.h>
#include <phat/representations/vector_list.h>
#include <phat/representations/sparse_pivot_column.h>
#include <phat/representations/heap_pivot_column.h>
#include <phat/representations/full_pivot_column.h>
#include <phat/representations/bit_tree_pivot_column.h>

// reduction algorithms
#include <phat/algorithms/twist_reduction.h>
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/spectral_sequence_reduction.h>

template<typename TPrecision, typename TPhatAlgorithm, typename TPhatDataStructure>
class ConeMapRipsPhat{


  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;



    typedef typename ConeMapFiltrationBuilder2<TPrecision>::NeighborMap NeighborMap;
    typedef typename ConeMapFiltrationBuilder2<TPrecision>::Neighbors Neighbors;
    typedef typename ConeMapFiltrationBuilder2<TPrecision>::NeighborsIterator NeighborsIterator;

    //typedef typename ConeMapFiltrationBuilder<TPrecision>::NeighborMap NeighborMap;
    //typedef typename ConeMapFiltrationBuilder<TPrecision>::Neighbors Neighbors;
    //typedef typename ConeMapFiltrationBuilder<TPrecision>::NeighborsIterator NeighborsIterator;



  private:
    GMRANeighborhood<TPrecision> &nh;
    MatrixXp dgm;
    MatrixXp death_centers;
    int filtrationSize;



  public:

    ConeMapRipsPhat( GMRANeighborhood<TPrecision> &nhood):nh(nhood), filtrationSize(0) {
    };



    void run( int maxD,
                  float radiusPercentile,
                  bool singleScale = false,
                  float radiusFactor = 1.01, bool store = false){

      if( radiusFactor <= 1){
        radiusFactor = 1.01;
      }
      //-- Instance to compute persistence homology at each scale  --//
      ConeMapHomology<TPrecision> homology;


      //-- Setup --//
      if(radiusPercentile < 0 ){
        radiusPercentile = 0;
      }

      //Find all leave nodes
      GMRATree<TPrecision> *gmra = nh.getTree();
      GMRAIDDecorator<TPrecision> ider;
      gmra->decorate(ider);
#ifdef SAVETOFILE
      int totalNumberOfNodes = ider.getMaxID();
#endif
      SetupRips setup;
      gmra->breadthFirstVisitor( &setup );

      NodeDistance<double> *dist =
         new CenterNodeDistance<double>( new EuclideanMetric<double>() );
      gmra->computeRadii(dist);
      delete dist;

      std::set< GMRAIDNode<TPrecision> * > &current = setup.leaves;


      TPrecision prevTime = 0;
      int scale = 0;
      int maxScale = 0;
      for(typename std::set<GMRAIDNode<TPrecision> *>::iterator it = current.begin();
                      it != current.end(); ++it){
          maxScale = std::max(maxScale, (*it)->getScale() );
      }


      //-- Run filtration at each scale --//

      //Birth simplicies to propagate across scale
      std::vector<IFiltration> fCones;
      IFiltration emptyCone;
      fCones.push_back( emptyCone );
      std::vector<Filtration> fScales;
      Filtration emptyF;
      double totalRipsTime = 0;
      double totalMapTime = 0;
      double totalNNTime = 0;
      while( current.size() > 1){


        clock_t t1 = clock();
#ifdef VERBOSE_DETAIL
        std::cout << std::endl << "Computing Rips at scale: " << scale << std::endl;
        std::cout << std::endl << "Computing Rips at scale: " << maxScale << std::endl;
        std::cout << "#Vertices : " << current.size() << std::endl;
#endif

        //-- Build filtration at current scale
        TPrecision radius = 0;
        if( singleScale ){
          radius = std::numeric_limits<TPrecision>::max() / 3;
        }else{
          radius = radiusFactor * getNthRadius(current, std::min( 1.0 +
                radiusPercentile * current.size(), current.size()-1.0) );
        }
#ifdef VERBOSE_DETAIL
        std::cout << "Radius Percentile: " << radiusPercentile << std::endl;
        std::cout << "Radius: " << radius << std::endl;
#endif

        for( typename std::set<GMRAIDNode<TPrecision> *>::iterator it = current.begin(); it
             != current.end(); ++it){
          (*it)->setStop( true );
        }

        //-- Compute neighborhood info
        Neighbors N;
        for(typename std::set<GMRAIDNode<TPrecision> *>::iterator it = current.begin(); it
            != current.end(); ++it){
          int index = (*it)->getID();
          NeighborMap &nl = N[index];

          typename GMRANeighborhood<TPrecision>::NeighborList nList;
          nh.neighbors(*it, 2*radius, nList);

          for(typename GMRANeighborhood<TPrecision>::NeighborListIterator nIt = nList.begin();
                          nIt != nList.end(); ++nIt){
              int nInd = dynamic_cast<GMRAIDNode<TPrecision> *>( nIt->second )->getID();
              nl[nInd] = nIt->first;
              N[nInd][index] = nIt->first;
          }
        }

        for(typename std::set<GMRAIDNode<TPrecision> *>::iterator it = current.begin();
                        it != current.end(); ++it){
          (*it)->setStop( false );
        }



        clock_t t2 =  clock();
        totalNNTime += ((double)t2 - t1) / CLOCKS_PER_SEC;
#ifdef VERBOSE_DETAIL
        std::cout << "Find neighbors time: " << ((double)t2 - t1) / CLOCKS_PER_SEC;
        std::cout << "N size : " << N.size() << std::endl;
        std::cout << std::endl;
        int nNeighbors = 0;
        for(NeighborsIterator it = N.begin(); it != N.end(); ++it){
          nNeighbors += it->second.size();
        }
        std::cout << "Average #neighbors: " << nNeighbors / (double) N.size();
        std::cout << std::endl;
#endif


        //-- Build filtration
        //Construct weighted filtration using the mapped filtration from the
        //previous scale
        IFiltration &lastCone = fCones.back();
        ConeMapFiltrationBuilder2<TPrecision> rips( lastCone);
        //ConeMapFiltrationBuilder<TPrecision> rips( );
        int maxSize = 100000000;
        fScales.push_back(emptyF);
        Filtration &currentF = fScales.back();
        rips.run( N, prevTime, 2*radius, maxSize, maxD, scale, currentF);

        clock_t t3 =  clock();
        totalRipsTime += ((double)t3 - t2) / CLOCKS_PER_SEC;
#ifdef VERBOSE_DETAIL
        std::cout << "Rips time: " << ((double)t3 - t2) / CLOCKS_PER_SEC << std::endl;
#endif



        //-- Compute homology at current scale
        TPrecision mapTime = radius;
        TPrecision maxTime = currentF.back().time;
#ifdef VERBOSE_DETAIL
        std::cout << "Rips max length: " << maxTime << std::endl;
        std::cout << std::endl << std::endl;
#endif




        //-- Construct map to next scale --//
        std::set< GMRAIDNode<TPrecision> * > next;
        if( !singleScale ){

          collectNodes( current, next, mapTime );
          if(next.size() == 1){
            break;
          }
          std::map<int, int> mapping;
          mapVertices( current, mapping, mapTime );

          //Construct mapping cone / cylinder
          fCones.push_back( emptyCone );
          IFiltration &cone = fCones.back();
          mapFiltration( currentF, mapping, cone, scale, maxTime, maxD );

          clock_t t4 =  clock();
          totalMapTime += ((double)t4 - t3) /CLOCKS_PER_SEC;
#ifdef VERBOSE_DETAIL
          std::cout << "Cone mapping filtration size: " << cone.size() << std::endl;
          std::cout << "Compute mapping time: " << ((double)t4 - t3) /CLOCKS_PER_SEC;
          std::cout << std::endl;
#endif
        }


        //-- Setup for next scale
        if(current.size() == next.size() ){
          radiusPercentile += 0.5 * (1-radiusPercentile);
          radiusFactor *= 1.5;
        }
        current = next;
        prevTime = maxTime;
        //prevEpsilon = epsilon;
        scale++;
        maxScale--;

        if(singleScale){
          break;
        }

      }


#ifdef VERBOSE
      clock_t t1 =  clock();
#endif

      //Create ids
      filtrationSize = 0;
      for(int i=0; i< fScales.size(); i++){
        filtrationSize += fScales[i].size();
      }
      for(int i=0; i< fCones.size(); i++){
        filtrationSize += fCones[i].size();
      }
      std::map<Simplex, int> ids;
      std::vector<double> times( filtrationSize );
      std::vector<double> scales( filtrationSize );
      std::vector<Simplex> filtration( filtrationSize );

#ifdef SAVETOFILE
      Eigen::VectorXi pointIDs(totalNumberOfNodes);
      for(int i=0; i<totalNumberOfNodes; i++){
        pointIDs[i] = -1;
      }
#endif

      for(int i=0; i< fScales.size(); i++){

         for( IFiltrationIterator it = fCones[i].begin(); it!=fCones[i].end(); ++it){
           const Simplex &s = it->first;
           int colIndex = ids.size();
           if(ids.find(s) == ids.end() ){
             ids[s] = colIndex;
             times[colIndex] = it->second;
             scales[colIndex] = i - 0.5;
             filtration[colIndex] = s;
#ifdef SAVETOFILE
             if(s.vertices.size() == 1){
               pointIDs( *s.vertices.begin() ) = colIndex;
             }
#endif
           }
         }

         for(FiltrationIterator it = fScales[i].begin(); it != fScales[i].end(); ++it){
           FiltrationEntry &fe = (*it);
           Simplex &s = fe.simplex;
           int colIndex = ids.size();
           if(ids.find(s) == ids.end() ){
             ids[s] = colIndex;
             times[colIndex] = fe.time;
             scales[colIndex] = i;
             filtration[colIndex] = s;
#ifdef SAVETOFILE
             if(s.vertices.size() == 1){
               pointIDs( *s.vertices.begin() ) = colIndex;
             }
#endif
           }
         }
      }

      filtrationSize = ids.size();

      //Build Phat boundary matrix
      phat::boundary_matrix< TPhatDataStructure > boundary_matrix;
      boundary_matrix.set_num_cols( filtrationSize );
      std::vector< phat::index > temp_col;
      for(int i=0; i<ids.size(); i++){
        const Simplex &s = filtration[i];
        int dim = s.vertices.size();
        boundary_matrix.set_dim(i, dim-1);
        if(dim == 1){
          boundary_matrix.set_col( i, temp_col);
        }
        else{
          getColumn(s, temp_col, ids);
          boundary_matrix.set_col( i, temp_col);
          temp_col.clear();
        }
      }

#ifdef VERBOSE
     clock_t t2 =  clock();
     std::cout << "---------------------------------------------" << std::endl;
     std::cout << "Filtration size: " <<  filtrationSize << std::endl;
     std::cout << "Ids size: " <<  ids.size() << std::endl;
     std::cout << "Construct Phat boundary matrix time: " << ((double)t2 - t1) /CLOCKS_PER_SEC;
     std::cout << std::endl;
#endif


#ifdef SAVETOFILE
      NodesToMatrix collectNodes( totalNumberOfNodes, gmra->getDataObject()->dimension() );
      gmra->breadthFirstVisitor( &collectNodes );
      EigenLinalg::LinalgIO<TPrecision>::writeMatrix(  "points.data", collectNodes.nodes);
      EigenLinalg::LinalgIO<int>::writeVector( "points-ids.data", pointIDs);
      Eigen::Map<VectorXp> v_scales( scales.data(), filtrationSize );
      EigenLinalg::LinalgIO<double>::writeVector( "scales.data", v_scales);
      boundary_matrix.save_ascii("boundary.txt");
#ifdef VERBOSE
     t2 =  clock();
#endif
#endif


      // define the object to hold the resulting persistence pairs
      phat::persistence_pairs pairs(times);

      // choose an algorithm (choice affects performance) and compute the persistence pair
      // (modifies boundary_matrix)
      phat::compute_persistence_pairs< TPhatAlgorithm >( pairs, boundary_matrix );

#ifdef VERBOSE
      clock_t t3 =  clock();
      std::cout << "Phat reduction time: " << ((double)t3 - t2) /CLOCKS_PER_SEC;
      std::cout << std::endl;
#endif
      //Store persistence diagram in matrix
      std::vector< phat::index > death_col;
      std::set< phat::index > death_pts;
      int index = 0;
      dgm = MatrixXp(4, pairs.get_num_pairs() );
      for(int i=0; i<pairs.get_num_pairs(); i++){
        int birth = pairs.get_pair( i ).first;
        int death = pairs.get_pair( i ).second;
        //if( times[death] - times[birth] > 0 ){
          dgm(0, index) = times[birth];
          dgm(1, index) = times[death];
          dgm(2, index) = scales[death];
          dgm(3, index) = boundary_matrix.get_dim( birth );
          boundary_matrix.get_col(death, death_col);
          for(int j=0; j<death_col.size(); j++){
            const Simplex &s = filtration[ death_col[j] ];
            death_pts.insert( s.vertices.begin(), s.vertices.end() );
          }
          ++index;
        //}
      }

      //dgm = dgm.leftCols(index);
      CollectNodes death_points(death_pts);
      gmra->breadthFirstVisitor( &death_points );
      
      int spatialDimension = gmra->getDataObject()->dimension();
      death_centers = MatrixXp(spatialDimension, pairs.get_num_pairs() );
      index = 0;
      for(int i=0; i<pairs.get_num_pairs(); i++){
        //int birth = pairs.get_pair( i ).first;
        int death = pairs.get_pair( i ).second;

        //if( times[death] - times[birth] > 0 ){
          boundary_matrix.get_col(death, death_col);
          VectorXp center = VectorXp::Zero(spatialDimension);
          int np = 0;
          for(int j=0; j<death_col.size(); j++){
            const Simplex &s = filtration[ death_col[j] ];
            for(std::set<int>::iterator sit = s.vertices.begin(); sit != s.vertices.end(); ++sit){
              //typename std::map< phat::index, VectorXp >::iterator it = death_points.points.find( *sit );
              //if( it != death_points.points.end() ){
              //  center += it->second;
              //}
              ++np;
              center += death_points.points[ *sit ];
            }
          }
          center /= np;
          death_centers.col(index) = center;
          ++index;
        //}
      }




#ifdef VERBOSE
      clock_t t4 =  clock();
      std::cout << "Store matrix time: " << ((double)t4 - t3) / CLOCKS_PER_SEC << std::endl;
      std::cout << "Total Rips time: " << totalRipsTime << std::endl;
      std::cout << "Total map time: " << totalMapTime << std::endl;
      std::cout << "Total NN time: " << totalNNTime << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
      std::cout << std::endl;
#endif


    };

    MatrixXp getDiagram(){
      return dgm;
    };

    int getFiltrationSize(){
        return filtrationSize;
    };

    MatrixXp getDeathCenters(){
      return death_centers;
    };



  private:

    class SetupRips : public Visitor<TPrecision>{

      public:

        std::set<GMRAIDNode<TPrecision> *> leaves;

        void visit(GMRANode<TPrecision> *node){
          //node->setNodeInfo( new RipsInfo() );
          if(node->getChildren().size() == 0){
            leaves.insert( dynamic_cast< GMRAIDNode<TPrecision> *>(node) );
          }
        };

    };

    class NodesToMatrix : public Visitor<TPrecision>{

      public:
        MatrixXp nodes;
        NodesToMatrix(int nNodes, int dim) : nodes(nNodes, dim) {
        };

        void visit(GMRANode<TPrecision> *node){
          int id = dynamic_cast< GMRAIDNode<TPrecision> *>(node)->getID();
          nodes.row( id ) = node->getCenter();
        };

    };

    class CollectNodes : public Visitor<TPrecision>{
      public:
        std::map< phat::index, VectorXp > points;
        std::set< phat::index > &nodes;
        CollectNodes(std::set< phat::index > &toCollect) : nodes(toCollect) {
        };

        void visit(GMRANode<TPrecision> *node){
          int id = dynamic_cast< GMRAIDNode<TPrecision> *>(node)->getID();
          if( nodes.find( id ) != nodes.end() ){
            points[id] = node->getCenter();
          }
        };

    };



    void getColumn( const Simplex &s, std::vector< phat::index > &col,
                    std::map<Simplex, int> &ids){
      std::list<Simplex> faces = s.getFaces();
      col.resize( s.vertices.size() );
      int index=0;
      for( std::list<Simplex>::iterator it = faces.begin(); it != faces.end();
           ++it, ++index){
        col[index] = ids[*it];
      }
      //col.sort();
      std::sort( col.begin(), col.end() );
    }



    void collectNodes( std::set< GMRAIDNode<TPrecision> * > &nodes,
                       std::set< GMRAIDNode<TPrecision> *>  &pnodes,
                       TPrecision r ){

      for(typename std::set< GMRAIDNode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++){
        pnodes.insert( getMappedNode(*it, r) );
      }

    };




    GMRAIDNode<TPrecision> *getMappedNode( GMRAIDNode<TPrecision> *node,
        TPrecision r){

      GMRAIDNode<TPrecision> *current = node;
      while( getRadius(current) <= r ){
      //if( getRadius(current) <= r ){
        GMRAIDNode<TPrecision> *parent = dynamic_cast<GMRAIDNode<TPrecision>*>( current->getParent() );
        if(parent != NULL){
          current = parent;
        }
        else{
          break;
        }
      }

      return current;
    };


    TPrecision getRadius(GMRANode<TPrecision> *node){
      return getParentRadius(node);


      TPrecision r = 0;
      GMRANode<TPrecision> *p = node->getParent();
      if(p == NULL){
        p = node;
      }

      std::vector< GMRANode<TPrecision>* > &kids = p->getChildren();

      NodeDistance<TPrecision> *dist = nh.getNodeDistance();
      for(unsigned int i=0; i<kids.size(); i++){
        r = std::max(r, dist->distance(kids[i], p) );
      }

      return r;
      //std::max( r, getRelativeRadius(node) );

    };




    TPrecision getParentRadius(GMRANode<TPrecision> *node){
      GMRANode<TPrecision> *p = node->getParent();
      if(p == NULL){
        return node->getRadius();
      }
      return p->getRadius();
    };


    TPrecision getMinRadius( std::set< GMRAIDNode<TPrecision> * > &nodes){
      TPrecision radius = std::numeric_limits<TPrecision>::max();
      for(typename std::set< GMRAIDNode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++){
        GMRANode<TPrecision> *n = *it;
        radius = std::min( getRadius(n), radius);
      }
      return radius;
    };

    TPrecision getNthRadius( std::set< GMRAIDNode<TPrecision> * > &nodes, int n){
      std::vector<double> radii(nodes.size());
      int index = 0;
      for(typename std::set< GMRAIDNode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++, index++){
        GMRANode<TPrecision> *n = *it;
        radii[index] =  getRadius(n);
      }
      std::sort( radii.begin(), radii.end() );
      return radii[n];
    };


    TPrecision getScaleMaxRadius( std::set< GMRAIDNode<TPrecision> * > &nodes, int scale){
      TPrecision radius = 0;
      for(typename std::set< GMRAIDNode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++ ){
        GMRANode<TPrecision> *n = *it;
        if(n->getScale() == scale){
          radius = std::max(radius, getRadius(n) );
        }
      }
      return radius;
    };


    void mapVertices( std::set< GMRAIDNode<TPrecision> *> &current,
                      std::map<int, int> &mapping, TPrecision r  ){

      for(typename std::set<GMRAIDNode<TPrecision> *>::iterator it = current.begin(); it !=
          current.end(); ++it){
        mapping[ (*it)->getID()] = getMappedNode( *it, r)->getID();
      }

    };


    //Map filtration according to mapping
    void mapFiltration( Filtration &filt, std::map<int, int> &mapping,
                        IFiltration &cone, int &scale, double &time,
                        int &maxDim ){


      for( FiltrationIterator it = filt.begin(); it != filt.end(); it++){

        Simplex &s = it->simplex;
        Simplex sNew;
        Simplex mapsTo;
        for(std::set<int>::iterator vIt = s.vertices.begin(); vIt !=
            s.vertices.end(); ++vIt){
          int v = *vIt;
          mapsTo.vertices.insert(v);
          if( v != mapping[v] ){
            sNew.vertices.insert( mapping[v] );
          }
        }

        if( sNew.vertices.size() > 0 ){
          cone[sNew] = time;
          if(   mapsTo.vertices.size() == 1  ){
            Simplex sCone = s;
            sCone.vertices.insert( *sNew.vertices.begin() );
            cone[sCone] = time;
          }
          else{
            joinSimplex( sNew, s, cone, scale, time, maxDim );
          }
        }
      }

    };


    void joinSimplex( Simplex &sMapped, Simplex &sOrig,
        IFiltration &cone, int &scale, double &time,
        int &maxDim){

      int nVertices = sMapped.vertices.size() + sOrig.vertices.size();
      bool keepGoing = true;
      if(nVertices <= maxDim + 2 ){
        Simplex sNew = sMapped;
        sNew.vertices.insert( sOrig.vertices.begin(), sOrig.vertices.end() );
        IFiltrationIterator fit = cone.find(sNew);
        keepGoing = fit == cone.end();
        if( keepGoing ){
          cone[sNew] = time;
        }
      }
      if( keepGoing){
        std::list<Simplex> f1 = sMapped.getFaces();
        for(std::list<Simplex>::iterator it = f1.begin(); it != f1.end(); ++it){
          joinSimplex( *it, sOrig, cone, scale, time, maxDim );
        }

        std::list<Simplex> f2 = sOrig.getFaces();
        for(std::list<Simplex>::iterator it = f2.begin(); it != f2.end(); ++it){
          joinSimplex( sMapped, *it, cone, scale, time, maxDim );
        }
      }

    };




};

#endif
