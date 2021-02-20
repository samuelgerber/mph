//author: Samuel Gerber


#ifndef MULTISCALEMAPPINGRIPS_H
#define MULTISCALEMAPPINGRIPS_H


#include "MappingFiltrationBuilder.h"
#include "GMRATree.h"
#include "GMRANeighborhood.h"
#include "EigenEuclideanMetric.h"

#include <map>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <vector>

#include <Eigen/Dense>

#include "RelativeGMRANeighborhood.h"
#include "MappingPersistentHomology.h"


template<typename TPrecision>
class MultiscaleMappingRips{


  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;


    typedef typename MappingPersistentHomology<TPrecision>::Event Event;
    typedef typename MappingPersistentHomology<TPrecision>::Events Events;
    typedef typename MappingPersistentHomology<TPrecision>::EventsIterator EventsIterator;

    typedef typename MappingPersistentHomology<TPrecision>::FiltrationEntryMap FiltrationEntryMap;
    typedef typename MappingPersistentHomology<TPrecision>::FiltrationEntryMapIterator FiltrationEntryMapIterator;

    typedef typename MappingFiltrationBuilder<TPrecision>::NeighborMap NeighborMap;
    typedef typename MappingFiltrationBuilder<TPrecision>::Neighbors Neighbors;




  private:


    GMRANeighborhood<TPrecision> &nh;
    //std::map< int, MappingFiltration > filtrations;
    std::map< int, MatrixXp > diagrams;

    int maximumFiltrationSize;
    int dim;



  public:

    MultiscaleMappingRips(GMRANeighborhood<TPrecision> &nhood) : nh(nhood){
      dim = nh.getTree()->getRoot()->getCenter().size();
    };



    void run( int maxD, float radiusPercentile, bool singleScale = false,
        float radiusFactor = 1.01 ){

      maximumFiltrationSize = 0;
      if( radiusFactor <= 1){
        radiusFactor = 1.01;
      }
      if(radiusPercentile < 0 ){
        radiusPercentile = 0;
      }



      //-- Setup --//

      //Instance to compute persistence homology at each scale//
      MappingPersistentHomology<TPrecision> homology;

      //Find all leave nodes
      SetupRips setup;
      GMRATree<TPrecision> *gmra = nh.getTree();
      gmra->breadthFirstVisitor( &setup );

      NodeDistance<double> *dist = new CenterNodeDistance<double>( new EuclideanMetric<double>() );
      gmra->computeRadii(dist);
      delete dist;

      std::set< GMRANode<TPrecision> * > current = setup.leaves;
      std::map< GMRANode<TPrecision> *, int > currentIndexes;
      //MatrixXp Xcurrent = getCenters(current, currentIndexes);
      getCenters(current, currentIndexes);

      std::vector<int> mapping;
      IMappingFiltration mapped;

      TPrecision prevTime = 0;
      //TPrecision prevEpsilon = 0;
      int scale = 0;
      int maxScale = 0;
      for(typename std::set<GMRANode<TPrecision> *>::iterator it = current.begin(); it
            != current.end(); ++it){
          maxScale = std::max(maxScale, (*it)->getScale() );
      }

      TPrecision maxRadius = gmra->getRoot()->getRadius();
      TPrecision currentRadius = 0;


      //-- Run filtration at each scale --//

      //Birth simplicies to propagate across scale
      std::set< GMRANode<TPrecision>* > activeNoMap;
      while( currentRadius < maxRadius){


#ifdef VERBOSE
        std::cout << std::endl << "Computing Rips at scale: " << scale << std::endl;
        std::cout << std::endl << "Computing Rips at scale: " << maxScale << std::endl;
        std::cout << "#Vertices : " << current.size() << std::endl;
        clock_t t1 = clock();
#endif

        //-- Build filtration at current scale

        //TPrecision radius = getNthRadius(current, std::min(current.size()/2.0,
        //     10.0) );
        //TPrecision radius = 1.2 * getMinRadius(current);
        TPrecision radius = 0;
        //TPrecision moveRadius = 1.5 * getMinRelativeRadius(current);
        //radius = std::max(moveRadius, radius);

        //TPrecision epsilon = radius;

        if( singleScale ){
          radius = std::numeric_limits<TPrecision>::max() / 3;
        }else{
          //radius = getScaleMaxRadius(current, maxScale);
          //radius = 1.5 * getMinRadius(current);
          radius = radiusFactor * getNthRadius(current, std::min( 1.0 +
                radiusPercentile * current.size(), current.size()-1.0) );
        }
#ifdef VERBOSE
        std::cout << "Radius Percentile: " << radiusPercentile << std::endl;
        std::cout << "Radius: " << radius << std::endl;
#endif
        currentRadius = radius;

        for(typename std::set<GMRANode<TPrecision> *>::iterator it = current.begin(); it
            != current.end(); ++it){
          (*it)->setStop( true );
        }

        //-- Compute neighborhood info
        Neighbors N(current.size());
        int index = 0;
        for(typename std::set<GMRANode<TPrecision> *>::iterator it = current.begin(); it
            != current.end(); ++it, ++index){

          NeighborMap &nl = N[index];

          typename GMRANeighborhood<TPrecision>::NeighborList nList;
          nh.neighbors(*it, 2*radius, nList);

          for(typename GMRANeighborhood<TPrecision>::NeighborListIterator nIt
              =nList.begin(); nIt != nList.end(); ++nIt){
              int nInd = currentIndexes[nIt->second];
              nl[nInd] = nIt->first;
              N[nInd][index] = nIt->first;
          }
        }

        for(typename std::set<GMRANode<TPrecision> *>::iterator it = current.begin(); it
            != current.end(); ++it){
          (*it)->setStop( false );
        }



#ifdef VERBOSE
        clock_t t2 =  clock();
        std::cout << "Find neighboors time: " << ((double)t2 - t1) / CLOCKS_PER_SEC;
        std::cout << std::endl;
#endif
        int nNeighbors = 0;
        for(unsigned int i=0; i<N.size(); i++){
          nNeighbors += N[i].size();
        }

#ifdef VERBOSE
        std::cout << "Average #neighbors: " << nNeighbors / (double) N.size();
        std::cout << std::endl;
#endif


        //-- Build filtration
        //Construct weighted filtration using the mapped filtration from the
        //previous scale
        MappingFiltrationBuilder<TPrecision> rips( current.size() );
        rips.setMappedMappingFiltration( mapped );
        int maxSize = 100000000;
        rips.run( N, prevTime, 2*radius, maxSize, maxD, scale);
        MappingFiltration filt = rips.getMappingFiltration();


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

        //Modified matrix reduction with mapping filtration
        if( maximumFiltrationSize < filt.size() ){
           maximumFiltrationSize = filt.size();
        };
        homology.run( filt, 2.0*mapTime, maxD );
        diagrams[scale] = homology.getDiagram();
        TPrecision maxTime = homology.getMaximalTime();

#ifdef VERBOSE
        clock_t t4 =  clock();
        std::cout << "Reduction time: " << ((double)t4 - t3) /CLOCKS_PER_SEC;
        std::cout << std::endl;
        std::cout << "MappingFiltration size: " << filt.size();
        std::cout << std::endl;
        std::cout << "Max Time: " << maxTime;
        std::cout << std::endl << std::endl;
#endif




        //-- Construct map to next scale --//
        std::set< GMRANode<TPrecision> * > next;
        std::map< GMRANode<TPrecision> *, int > nextIndexes;

        //Compute vertex mapping
        if( !singleScale ){

          //Make sure birth simplicies ar kept intact
          std::set< GMRANode<TPrecision> *> noMap;
          std::vector< GMRANode<TPrecision> * > currentNodes( currentIndexes.size() );
          for(typename std::map<GMRANode<TPrecision> *, int>::iterator it =
              currentIndexes.begin(); it!=currentIndexes.end(); ++it){
            currentNodes[it->second] = it->first;
            GMRANode<TPrecision> *p = it->first->getParent();
            if( p != NULL){
              p = p->getParent();
              while( p !=NULL ){
                noMap.insert(p);
                p = p->getParent();
              }
            }
          }

          FiltrationEntryMap &active = homology.getActive();
          for(FiltrationEntryMapIterator it = active.begin(); it !=active.end(); ++it){
            std::set< GMRANode<TPrecision>* > parents;
            MappingFiltrationEntry &entry = it->second;
            MappingSimplex &simplex = entry.simplex;
            std::cout << "Mapping birth time: " << entry.time << std::endl;
            std::cout << "Mapping simplex scale: " << simplex.getScale() << std::endl;
            for( std::set<int>::iterator it = simplex.vertices.begin(); it
                   !=simplex.vertices.end(); ++it){
              GMRANode<TPrecision> *p = currentNodes[*it]->getParent();
              if( p != NULL ){
                typename std::set< GMRANode<TPrecision>* >::iterator pIt = parents.find(p);
                activeNoMap.insert( p );
                /*
                if( pIt == parents.end() ){
                  parents.insert( p );
                }
                else{
                  activeNoMap.insert( p );
                }
                */
              }
            }
          }
          std::cout << "#ActiveNoMap: " << activeNoMap.size() << std::endl;
          noMap.insert( activeNoMap.begin(), activeNoMap.end() );


          collectNodes( current, next, mapTime, noMap );
          //make sure active birth simplicies are mapped to individual simnplicies


          //Xnext = getCenters( next, nextIndexes );
          getCenters( next, nextIndexes );
          mapVertices( currentIndexes, nextIndexes, mapping, mapTime, noMap );

          //make sure active birth simplicies are preserved across the mapping
          /*
          std::vector< GMRANode<TPrecision> * > nextNodes( nextIndexes.size() );
          for(typename std::map< GMRANode<TPrecision> *, int>::iterator it =
              nextIndexes.begin(); it!=nextIndexes.end(); ++it){
            nextNodes[it->second] = it->first;
          }
          std::vector< GMRANode<TPrecision> * > currentNodes( currentIndexes.size() );
          for(typename std::map<GMRANode<TPrecision> *, int>::iterator it =
              currentIndexes.begin(); it!=currentIndexes.end(); ++it){
            currentNodes[it->second] = it->first;
          }

          FiltrationEntryMap &active = homology.getActive();
          for( FiltrationEntryMapIterator it = active.begin(); it !=active.end(); ++it){
             MappingSimplex &simplex = it->second;
             std::set<int> mapsTo;
             for( std::set<int>::iterator it = simplex.vertices.begin(); it
                 !=simplex.vertices.end(); ++it){
               int mId = mapping[ *it ];
               std::set<int>::iterator mIt = mapsTo.find( mId );
               if( mIt == mapsTo.end() ){
                 mapsTo.insert(mId);
               }
               else{
                 std::cout << "Moving Node" << std::endl;
                 //Insert new GMRA Node
                 GMRANode<TPrecision> *n = nextNodes[ mId ];
                 GMRANode<TPrecision> *np = n->getParent();
                 GMRANode<TPrecision> *c = currentNodes[ *it ];
                 GMRANode<TPrecision> *cp = c->getParent();
                 if( cp != NULL){
                   cp->removeChild(c);
                 }
                 else{
                   std::cout << "nsfw1" << std::endl;
                 }
                 if( np != NULL ){
                   //hacky destroys internal structure of specialized GMRA
                   //trees (e.g. IPCATree)
                   np->getChildren().push_back( c );
                   c->setParent( np );
                   next.insert( c );
                   int index = nextIndexes.size();
                   nextIndexes[ c ] = index;
                   mapping[ *it ] = index;
                   mapsTo.insert(index);
                   nextNodes.push_back( c );
                 }
                 else{
                   std::cout << "nsfw2" << std::endl;
                 }

               }

             }

          }
          */


        }
        else{
          mapTime = std::numeric_limits<TPrecision>::max()/3;
        }


#ifdef VERBOSE
        clock_t t5 =  clock();
        std::cout << "Compute mapping time: " << ((double)t5 - t4) /CLOCKS_PER_SEC;
        std::cout << std::endl;
#endif


        //Map filtration to next scale scale
        if( !singleScale ){
          mapped.clear();
          FiltrationEntryMap &active = homology.getActive();
          mapMappingFiltration(filt, mapping, mapped, scale, active);
        }

#ifdef VERBOSE
        clock_t t6 =  clock();
        std::cout << "Map filtration  time: " << ((double)t6 - t5) /CLOCKS_PER_SEC;
        std::cout << std::endl;
#endif


        //-- Setup for next scale
        if(current.size() == next.size() ){
          radiusPercentile += 0.5 * (1-radiusPercentile);
          radiusFactor *= 1.5;
        }
        current = next;
        currentIndexes = nextIndexes;
        prevTime = maxTime;
        //prevEpsilon = epsilon;
        scale++;
        maxScale--;

        if(singleScale){
          break;
        }

      }

    };

    int getMaximumFiltrationSize(){
      return maximumFiltrationSize;
    };



    std::map<int, MatrixXp > &getDiagrams(){
      return diagrams;
    };



  private:

    class SetupRips : public Visitor<TPrecision>{

      public:

        std::set<GMRANode<TPrecision> *> leaves;

        void visit(GMRANode<TPrecision> *node){
          //node->setNodeInfo( new RipsInfo() );
          if(node->getChildren().size() == 0){
            leaves.insert(node);
          }
        };

    };



    void collectNodes(std::set< GMRANode<TPrecision> * > &nodes,
        std::set<GMRANode<TPrecision> *> &pnodes, TPrecision r,
        std::set< GMRANode<TPrecision>* > &noMap){

      for(typename std::set< GMRANode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++){
        pnodes.insert( getMappedNode(*it, r, noMap) );
      }

    };




    GMRANode<TPrecision> *getMappedNode( GMRANode<TPrecision> *node,
        TPrecision r, std::set< GMRANode<TPrecision>* > &noMap){

      GMRANode<TPrecision> *current = node;
      //while( getRadius(current) <= r ){
      if( getRadius(current) <= r ){
        GMRANode<TPrecision> *parent = current->getParent();
        if(parent != NULL){
          if( noMap.find(parent) == noMap.end() ){
            current = parent;
          }
        }
        //else{
        //  break;
        //}
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

/*
 *

    TPrecision getRelativeRadius(GMRANode<TPrecision> *from,
        GMRANode<TPrecision> *to){

      GMRANode<TPrecision> *p = to->getParent();
      if(p == NULL){
        return std::numeric_limits<TPrecision>::max();
      }

      static EuclideanMetric<TPrecision> l2norm;
      VectorXp x1 = p->getCenter();
      VectorXp x2 = from->getCenter();
      TPrecision r = l2norm.distance(x1, x2);

      return r;
    };


    TPrecision getRelativeRadius(GMRANode<TPrecision> *node){

      GMRANode<TPrecision> *p = node->getParent();
      if(p == NULL){
        return std::numeric_limits<TPrecision>::max();
      }

      TPrecision r = 0;

      std::vector< GMRANode<TPrecision>* > &kids = p->getChildren();

      static EuclideanMetric<TPrecision> l2norm;
      VectorXp x1 = p->getCenter();
      for(int i=0; i<kids.size(); i++){
          VectorXp x2 = kids[i]->getCenter();
          TPrecision tmp = l2norm.distance(x1, x2);
          r = std::max(r, tmp);
      }

      return r;
    };




    TPrecision getMinRadius( std::set< GMRANode<TPrecision> * > &nodes){
      TPrecision r = std::numeric_limits<TPrecision>::max();
      for(typename std::set< GMRANode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++){
        GMRANode<TPrecision> *n = *it;
        r = std::min(r, getRadius(n) );
      }
      return r;
    };



    TPrecision getMinRelativeRadius( std::set< GMRANode<TPrecision> * > &nodes){
      TPrecision r = std::numeric_limits<TPrecision>::max();
      for(typename std::set< GMRANode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++){
        GMRANode<TPrecision> *n = *it;
        r = std::min(r, getRelativeRadius(n) );
      }
      return r;
    };



*/


    TPrecision getMinRadius( std::set< GMRANode<TPrecision> * > &nodes){
      TPrecision radius = std::numeric_limits<TPrecision>::max();
      for(typename std::set< GMRANode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++){
        GMRANode<TPrecision> *n = *it;
        radius = std::min(getRadius(n), radius);
      }
      return radius;
    };

    TPrecision getNthRadius( std::set< GMRANode<TPrecision> * > &nodes, int n){
      std::vector<double> radii(nodes.size());
      int index = 0;
      for(typename std::set< GMRANode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++, index++){
        GMRANode<TPrecision> *n = *it;
        radii[index] =  getRadius(n);
      }
      std::sort( radii.begin(), radii.end() );
      return radii[n];
    };


    TPrecision getScaleMaxRadius( std::set< GMRANode<TPrecision> * > &nodes, int scale){
      TPrecision radius = 0;
      for(typename std::set< GMRANode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++ ){
        GMRANode<TPrecision> *n = *it;
        if(n->getScale() == scale){
          radius = std::max(radius, getRadius(n) );
        }
      }
      return radius;
    };





    void mapVertices( std::map< GMRANode<TPrecision> *, int > &X,
                      std::map< GMRANode<TPrecision> *, int > &Y,
                      std::vector<int> &mapping, TPrecision r,
                      std::set< GMRANode<TPrecision>* > &noMap){

      mapping.clear();
      mapping.resize(X.size());

      for(typename std::map<GMRANode<TPrecision> *, int>::iterator it = X.begin(); it !=
          X.end(); ++it){
        mapping[it->second] = Y[ getMappedNode(it->first, r, noMap) ];
      }

    };


    void getCenters(std::set<GMRANode<TPrecision> *> current,
        std::map<GMRANode<TPrecision> *, int> &indexes){
      //Collect centers for this scale
      //MatrixXp Xs(dim, current.size());
      int index= 0;
      for(typename std::set<GMRANode<TPrecision> *>::iterator it = current.begin(); it !=
          current.end(); ++it, ++index){
        GMRANode<TPrecision> *node = (*it);
        indexes[node] = index;
        //Xs.col(index) = node->getCenter();
      }
      //return Xs;
    };

/*
    MatrixXp getCenters(std::set<GMRANode<TPrecision> *> current,
        std::map<GMRANode<TPrecision> *, int> &indexes){
      using namespace FortranLinalg;
      //Collect centers for this scale
      MatrixXp Xs(dim, current.size());
      int index= 0;
      for(typename std::set<GMRANode<TPrecision> *>::iterator it = current.begin(); it !=
          current.end(); ++it, ++index){
        GMRANode<TPrecision> *node = (*it);
        indexes[node] = index;
        Xs.col(index) = node->getCenter();
      }
      return Xs;
    };

*/

    //Map filtration according to mapping
    void mapMappingFiltration(MappingFiltration &filt, std::vector<int> &mapping,
              IMappingFiltration &mapped, int scale, FiltrationEntryMap &active){

      for( MappingFiltrationIterator it = filt.begin(); it != filt.end(); it++){
        MappingSimplex sNew(scale);

        //
        MappingSimplex &s = it->simplex;

        for(std::set<int>::iterator vIt = s.vertices.begin(); vIt !=
            s.vertices.end(); ++vIt){
          int v = *vIt;
          sNew.vertices.insert( mapping[v] );
        }
        //if(s.vertices.size() == sNew.vertices.size() ){
        //  s.setMapsTo(sNew.vertices);
        //}

        if(sNew.vertices.size() > 0 ){

          IMappingFiltrationIterator it2 = mapped.find(sNew);
          if(it2 == mapped.end()){
            sNew.setScale(s.getScale());
            mapped[sNew] = it->time;
          }
          else{
            const MappingSimplex &s2 = it2->first;
            if(s.getScale() < s2.getScale()){
              mapped.erase(sNew);
              sNew.setScale( s.getScale() );
              mapped[sNew] = it->time;
            }
            else if(s.getScale() == s2.getScale() ){
              it2->second =  std::min(it->time, it2->second);
            }
            //else leave old one
          }
        }
        else{
#ifdef VERBOSE
          std::cout << "Huuuummmmm" << std::endl;
#endif
        }
      }



      //Overwrite with active birth simplicies
      for( FiltrationEntryMapIterator it = active.begin(); it != active.end(); ++it){
        MappingFiltrationEntry &entry = it->second;
        MappingSimplex &s = entry.simplex;
        MappingSimplex sNew( s.getScale() );
        for(std::set<int>::iterator vIt = s.vertices.begin(); vIt !=
            s.vertices.end(); ++vIt){
          int v = *vIt;
          sNew.vertices.insert( mapping[v] );
        }
        if( sNew.vertices.size()  != s.vertices.size() ){
          std::cout << "--------- Active map error" << std::endl;
        }
        IMappingFiltrationIterator it2 = mapped.find(sNew);
        if(it2 == mapped.end()){
          std::cout << "--------- Active map error" << std::endl;
          mapped[sNew] = entry.time;
        }
        else{
          std::cout << "Propagated mapping birth time: " << it2->second << std::endl;;
          std::cout << "Propagated mapping birth scale: " << it2->first.getScale() << std::endl;;
          it2->second =  entry.time;
          std::cout << "Corrected mapping birth time: " << entry.time << std::endl;;
          std::cout << "Corrected mapping birth scale: " << s.getScale() << std::endl;;
        }
      }


      //Make sure the order of simplicies is preserved. It can happen that a
      //lower order MappingSimplex has a larger time stamp
      makeConsistent(mapped);
#ifdef VERBOSE
      std::cout << "Mapped size: " << mapped.size() << std::endl;
#endif

    };



    void makeConsistent(IMappingFiltration &ifilt){

      for(IMappingFiltrationIterator it = ifilt.begin(); it != ifilt.end(); ++it){

        MappingSimplex s = it->first;
        TPrecision tMax = it->second;

        if(s.vertices.size() > 1){
          //check maximal time of face addition
          std::list<MappingSimplex> sList = s.getFaces();
          for(std::list<MappingSimplex>::iterator sIt = sList.begin(); sIt !=
              sList.end(); ++sIt){

            MappingSimplex &f = *sIt;

            IMappingFiltrationIterator fIt = ifilt.find(f);
            if(fIt != ifilt.end() ){
              const MappingSimplex &f2 = fIt->first;
              if( f2.getScale() <= s.getScale() ){
                if(tMax < fIt->second){
                  tMax = fIt->second;
                }
              }
              else{
#ifdef VERBOSE
                std::cout << "-------- woof" << std::endl;
#endif
              }
            }
            else{
#ifdef VERBOSE
              std::cout << "---------- bark" << std::endl;
#endif
            }
          }

          it->second = tMax;
        }
      }

    };




/*
    //----Debug ----//

    struct Triple{
      Triple(int s, int index1, int index2):scale(s), i1(index1), i2(index2){};
      int scale;
      int i1;
      int i2;
    };


    //Debug method
    void writeEdges(IMappingFiltration &ifilt, std::string filename){
      using namespace FortranLinalg;
      std::list< Triple  > edges;
      for(IMappingFiltrationIterator it = ifilt.begin(); it != ifilt.end(); ++it){
        const MappingSimplex &s = it->first;
        if(s.vertices.size() == 2){
          int i1 = *s.vertices.begin();
          int i2 = *s.vertices.rbegin();
          Triple edge(s.getScale(), i1, i2);
          edges.push_back(edge);
        }
      }

      MatrixXp E(3, edges.size());
      int index = 0;
      for(typename std::list< Triple >::iterator it = edges.begin(); it!=
          edges.end(); ++it, ++index){
        E(0, index) = it->i1;
        E(1, index) = it->i2;
        E(2, index) = it->scale;
      }
    };


    //Debug method
    void writeEdges(MappingFiltration &filt, std::string filename){
      using namespace FortranLinalg;
      std::list< Triple  > edges;
      for(MappingFiltrationIterator it = filt.begin(); it != filt.end(); ++it){
        const MappingSimplex &s = it->simplex;
        if(s.vertices.size() == 2){
          int i1 = *s.vertices.begin();
          int i2 = *s.vertices.rbegin();
          Triple edge(s.getScale(), i1, i2);
          edges.push_back(edge);
        }
      }

      MatrixXp E(3, edges.size());
      int index = 0;
      for(typename std::list< Triple >::iterator it = edges.begin(); it!=
          edges.end(); ++it, ++index){
        E(0, index) = it->i1;
        E(1, index) = it->i2;
        E(2, index) = it->scale;
      }
    };
*/


};

#endif
