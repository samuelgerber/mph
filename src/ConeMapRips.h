//author: Samuel Gerber


#ifndef CONEMAPRIPS_H
#define CONEMAPRIPS_H


#include "GMRATree.h"
#include "GMRANeighborhood.h"
#include "RelativeGMRANeighborhood.h"
#include "GMRAIDNode.h"

#include "ConeMapFiltrationBuilder.h"
#include "ConeMapHomology.h"

#include <map>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <vector>

template<typename TPrecision>
class ConeMapRips{


  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;


    typedef typename ConeMapHomology<TPrecision>::Event Event;
    typedef typename ConeMapHomology<TPrecision>::Events Events;
    typedef typename ConeMapHomology<TPrecision>::EventsIterator EventsIterator;

    typedef typename ConeMapFiltrationBuilder<TPrecision>::NeighborMap NeighborMap;
    typedef typename ConeMapFiltrationBuilder<TPrecision>::Neighbors Neighbors;
    typedef typename ConeMapFiltrationBuilder<TPrecision>::NeighborsIterator NeighborsIterator;




  private:
    GMRANeighborhood<TPrecision> &nh;
    MatrixXp dgm;
    int filtrationSize;



  public:

    ConeMapRips( GMRANeighborhood<TPrecision> &nhood):nh(nhood), filtrationSize(0) {
    };



    void run( int maxD,
                  float radiusPercentile,
                  bool singleScale = false,
                  float radiusFactor = 1.01 ){

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
      while( current.size() > 1){


#ifdef VERBOSE
        std::cout << std::endl << "Computing Rips at scale: " << scale << std::endl;
        std::cout << std::endl << "Computing Rips at scale: " << maxScale << std::endl;
        std::cout << "#Vertices : " << current.size() << std::endl;
        clock_t t1 = clock();
#endif

        //-- Build filtration at current scale
        TPrecision radius = 0;
        if( singleScale ){
          radius = std::numeric_limits<TPrecision>::max() / 3;
        }else{
          radius = radiusFactor * getNthRadius(current, std::min( 1.0 +
                radiusPercentile * current.size(), current.size()-1.0) );
        }
#ifdef VERBOSE
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



#ifdef VERBOSE
        clock_t t2 =  clock();
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
        ConeMapFiltrationBuilder<TPrecision> rips;
        int maxSize = 100000000;
        rips.run( N, prevTime, 2*radius, maxSize, maxD, scale);
        Filtration filt = rips.getFiltration();

#ifdef VERBOSE
        clock_t t3 =  clock();
        std::cout << "Rips time: " << ((double)t3 - t2) / CLOCKS_PER_SEC << std::endl;
#endif



        //-- Compute homology at current scale
        TPrecision mapTime = radius;
        if( rips.isTruncated() ){
          TPrecision ripsTime = rips.getMaxTime();
#ifdef VERBOSE
          std::cout << "Rips max length: " << ripsTime << std::endl;
#endif
          mapTime = std::min( ripsTime/2, radius );
        }

        homology.run( filt, maxD, scale );
        TPrecision maxTime = homology.getMaximalTime();

#ifdef VERBOSE
        clock_t t4 =  clock();
        std::cout << "Reduction time: " << ((double)t4 - t3) /CLOCKS_PER_SEC;
        std::cout << std::endl;
        std::cout << "Filtration size: " << filt.size();
        std::cout << std::endl;
        std::cout << "Max Time: " << maxTime;
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
          IFiltration cone;
          mapFiltration( filt, mapping, cone, scale, maxTime, maxD );

#ifdef VERBOSE
          clock_t t5 =  clock();
          std::cout << "Cone mapping filtration size: " << cone.size() << std::endl;
          std::cout << "Compute mapping time: " << ((double)t5 - t4) /CLOCKS_PER_SEC;
          std::cout << std::endl;
#endif

          homology.run( cone, maxD, scale );

#ifdef VERBOSE
          clock_t t6 =  clock();
          std::cout << "Cone map homology reduction  time: ";
          std::cout << ((double)t6 - t5) / CLOCKS_PER_SEC;
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

      dgm = homology.getDiagram();
      filtrationSize = homology.getFiltrationSize();

    };



    MatrixXp getDiagram(){
      return dgm;
    };



    int getFiltrationSize(){
        return filtrationSize;
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







    void collectNodes( std::set< GMRAIDNode<TPrecision> * > &nodes,
                       std::set< GMRAIDNode<TPrecision> *> &pnodes,
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
          //sNew.setMapped( true );
          cone[sNew] = time;
          //long nVertices = s.vertices.size();
          //if( nVertices < maxDim + 2 ){
            //Add cone simplicies
            if(   mapsTo.vertices.size() == 1  ){
              Simplex sCone = s;
              //sCone.setMapped(true);
              sCone.vertices.insert( *sNew.vertices.begin() );
              cone[sCone] = time;
            }
            else{
              addCones( sNew, s, cone, scale, time, maxDim );
            }
          //}
          //if(mapsTo.vertices.size() != sNew.vertices.size() ){
            //mapsTo.setMapped( true );
            //cone[mapsTo] = time;
          //}
        }
      }

    };

    void addCones( Simplex &sMapped, Simplex &sOrig,
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
          //sNew.setMapped( true );
          cone[sNew] = time;
        }
      }
      if( keepGoing){
        std::list<Simplex> f1 = sMapped.getFaces();
        for(std::list<Simplex>::iterator it = f1.begin(); it != f1.end(); ++it){
          addCones( *it, sOrig, cone, scale, time, maxDim );
        }

        std::list<Simplex> f2 = sOrig.getFaces();
        for(std::list<Simplex>::iterator it = f2.begin(); it != f2.end(); ++it){
          addCones( sMapped, *it, cone, scale, time, maxDim );
        }
      }

    };




};

#endif
