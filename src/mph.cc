#ifndef NULL
#define NULL 0
#endif

#define R_NO_REMAP

#define VERBOSE
//#define SAVETOFILE

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdio.h>


#include "GMRATree.h"
#include "GMRANeighborhood.h"
#include "EigenL1Metric.h"
#include "EigenSquaredEuclideanMetric.h"
#include "EigenEuclideanMetric.h"
#include "EuclideanMetric.h"
#include "WassersteinNodeDistance.h"

#include "Distance.h"
#include "DenseMatrix.h"

#include "MultiscaleMappingRips.h"
#include "MultiscaleMappingRips2.h"
#include "ConeMapRips.h"
#include "ConeMapRipsPhat.h"

#include "FiltrationBuilder.h"
#include "PersistentHomology.h"
#include "QuadPersistentHomology.h"
#include "Grid3dNeighbors.h"

#include "NodeDistance.h"




extern "C" {
  //Utility functions
  NodeDistance<double> *getNodeDistance(int distType){
    NodeDistance<double> *dist = NULL;
    if(distType == 1){
      dist = new CenterNodeDistance<double>( new EuclideanMetric<double>() );
    }
    else if(distType == 2){
      dist = new CenterNodeDistance<double>( new L1Metric<double>() );
    }
    else if(distType == 3){
      dist = new CenterNodeDistance<double>( new SquaredEuclideanMetric<double>() );
    }
    else if(distType == 4){
      dist = new WassersteinNodeDistance<double>();
    }
    return dist;
  };


  /**/
  /* Multiscale Topology Calls */

  //Multiscale Rips from GMRA tree
  SEXP runMultiscaleRips( GMRANeighborhood<double> &nh,  int maxD, float percentile,
                          bool single, float radiusFactor, int type, int algType, int dsType){
    using namespace Eigen;

    MatrixXd dgm;
    MatrixXd centers(0 , 0);
    int fSize;
    if(type == 0){
      ConeMapRips<double> mrips(nh);
      mrips.run(maxD, percentile, single, radiusFactor);
      dgm = mrips.getDiagram();
      fSize = mrips.getFiltrationSize();
    }
    else{
      switch(dsType){
        case 1:
          switch(algType){
            case 1:
              {ConeMapRipsPhat<double, phat::twist_reduction, phat::vector_vector> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 2:
              {ConeMapRipsPhat<double, phat::standard_reduction, phat::vector_vector> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 3:
              {ConeMapRipsPhat<double, phat::row_reduction, phat::vector_vector> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 4:
              {ConeMapRipsPhat<double, phat::chunk_reduction, phat::vector_vector> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 5:
              {ConeMapRipsPhat<double, phat::spectral_sequence_reduction, phat::vector_vector> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
          }
          break;
        case 2:
          switch(algType){
            case 1:
              {ConeMapRipsPhat<double, phat::twist_reduction, phat::vector_heap> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 2:
              {ConeMapRipsPhat<double, phat::standard_reduction, phat::vector_heap> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 3:
              {ConeMapRipsPhat<double, phat::row_reduction, phat::vector_heap> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 4:
              {ConeMapRipsPhat<double, phat::chunk_reduction, phat::vector_heap> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 5:
              {ConeMapRipsPhat<double, phat::spectral_sequence_reduction, phat::vector_heap> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
          }
          break;
        case 3:
          switch(algType){
            case 1:
              {ConeMapRipsPhat<double, phat::twist_reduction, phat::vector_set> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 2:
              {ConeMapRipsPhat<double, phat::standard_reduction, phat::vector_set> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 3:
              {ConeMapRipsPhat<double, phat::row_reduction, phat::vector_set> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 4:
              {ConeMapRipsPhat<double, phat::chunk_reduction, phat::vector_set> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 5:
              {ConeMapRipsPhat<double, phat::spectral_sequence_reduction, phat::vector_set> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
          }
          break;
        case 4:
          switch(algType){
            case 1:
              {ConeMapRipsPhat<double, phat::twist_reduction, phat::vector_list> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 2:
              {ConeMapRipsPhat<double, phat::standard_reduction, phat::vector_list> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              fSize = mrips.getFiltrationSize();}
              break;
            case 3:
              {ConeMapRipsPhat<double, phat::row_reduction, phat::vector_list> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 4:
              {ConeMapRipsPhat<double, phat::chunk_reduction, phat::vector_list> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 5:
              {ConeMapRipsPhat<double, phat::spectral_sequence_reduction, phat::vector_list> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
          }
          break;
        case 5:
          switch(algType){
            case 1:
              {ConeMapRipsPhat<double, phat::twist_reduction, phat::heap_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 2:
              {ConeMapRipsPhat<double, phat::standard_reduction, phat::heap_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 3:
              {ConeMapRipsPhat<double, phat::row_reduction, phat::heap_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 4:
              {ConeMapRipsPhat<double, phat::chunk_reduction, phat::heap_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 5:
              {ConeMapRipsPhat<double, phat::spectral_sequence_reduction, phat::heap_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
          }
          break;
        case 6:
          switch(algType){
            case 1:
              {ConeMapRipsPhat<double, phat::twist_reduction, phat::full_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 2:
              {ConeMapRipsPhat<double, phat::standard_reduction, phat::full_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 3:
              {ConeMapRipsPhat<double, phat::row_reduction, phat::full_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 4:
              {ConeMapRipsPhat<double, phat::chunk_reduction, phat::full_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 5:
              {ConeMapRipsPhat<double, phat::spectral_sequence_reduction, phat::full_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
          }
          break;
        case 7:
          switch(algType){
            case 1:
              {ConeMapRipsPhat<double, phat::twist_reduction, phat::sparse_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 2:
              {ConeMapRipsPhat<double, phat::standard_reduction, phat::sparse_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 3:
              {ConeMapRipsPhat<double, phat::row_reduction, phat::sparse_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 4:
              {ConeMapRipsPhat<double, phat::chunk_reduction, phat::sparse_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 5:
              {ConeMapRipsPhat<double, phat::spectral_sequence_reduction, phat::sparse_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
          }
          break;
        case 8:
           switch(algType){
            case 1:
              {ConeMapRipsPhat<double, phat::twist_reduction, phat::bit_tree_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 2:
              {ConeMapRipsPhat<double, phat::standard_reduction, phat::bit_tree_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 3:
              {ConeMapRipsPhat<double, phat::row_reduction, phat::bit_tree_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 4:
              {ConeMapRipsPhat<double, phat::chunk_reduction, phat::bit_tree_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
            case 5:
              {ConeMapRipsPhat<double, phat::spectral_sequence_reduction, phat::bit_tree_pivot_column> mrips(nh);
              mrips.run(maxD, percentile, single, radiusFactor);
              dgm = mrips.getDiagram();
              centers = mrips.getDeathCenters();
              fSize = mrips.getFiltrationSize();}
              break;
          }
          break;
      }
    }


    SEXP list;
    PROTECT( list = Rf_allocVector(VECSXP, 3) );

    SEXP Rdgm;
    PROTECT( Rdgm = Rf_allocMatrix(REALSXP, dgm.rows(), dgm.cols()) );
    memcpy( REAL(Rdgm), dgm.data(), dgm.cols()*dgm.rows()*sizeof(double) );
    SET_VECTOR_ELT(list, 0, Rdgm);

    SEXP RfSize;
    PROTECT(RfSize = Rf_allocVector(INTSXP, 1));
    INTEGER(RfSize)[0] = fSize;
    SET_VECTOR_ELT(list, 1, RfSize);

    SEXP Rcenters;
    PROTECT( Rcenters = Rf_allocMatrix(REALSXP, centers.rows(), centers.cols()) );
    memcpy( REAL(Rcenters), centers.data(), centers.cols()*centers.rows()*sizeof(double) );
    SET_VECTOR_ELT(list, 2, Rcenters);

    UNPROTECT(4);

    return list;

  };


  SEXP runMultiscaleRipsRisky1( GMRANeighborhood<double> &nh,  int maxD, float percentile,
                               bool single, float radiusFactor){
    using namespace Eigen;

    MultiscaleMappingRips<double> mrips(nh);
    mrips.run(maxD, percentile, single, radiusFactor);
    std::map<int, MatrixXd> &dgms = mrips.getDiagrams();

    int size = 0;
    for(std::map<int, MatrixXd >::iterator it = dgms.begin(); it !=
        dgms.end(); ++it){
      size += it->second.cols();
    };

    MatrixXd dgm(4, size);
    int index = 0;
    for(std::map<int, MatrixXd >::iterator it = dgms.begin(); it !=
        dgms.end(); ++it){
      MatrixXd &dg = it->second;
      for(unsigned int i=0; i< dg.cols(); i++, index++){
        dgm(0, index) = dg(0, i);
        dgm(1, index) = dg(1, i);
        dgm(2, index) = it->first;
        dgm(3, index) = dg(2, i);
      }
    }
    int fSize = mrips.getMaximumFiltrationSize();

    SEXP list;
    PROTECT( list = Rf_allocVector(VECSXP, 2) );

    SEXP Rdgm;
    PROTECT( Rdgm = Rf_allocMatrix(REALSXP, dgm.rows(), dgm.cols()) );
    memcpy( REAL(Rdgm), dgm.data(), dgm.cols()*dgm.rows()*sizeof(double) );
    SET_VECTOR_ELT(list, 0, Rdgm);

    SEXP RfSize;
    PROTECT(RfSize = Rf_allocVector(INTSXP, 1));
    INTEGER(RfSize)[0] = fSize;
    SET_VECTOR_ELT(list, 1, RfSize);

    UNPROTECT(3);

    return list;

  };

  SEXP runMultiscaleRipsRisky2( GMRANeighborhood<double> &nh,  int maxD, float percentile,
                               bool single, float radiusFactor){
    using namespace Eigen;

    MultiscaleMappingRips2<double> mrips(nh);
    mrips.run(maxD, percentile, single, radiusFactor);
    std::map<int, MatrixXd> &dgms = mrips.getDiagrams();

    int size = 0;
    for(std::map<int, MatrixXd >::iterator it = dgms.begin(); it !=
        dgms.end(); ++it){
      size += it->second.cols();
    };

    MatrixXd dgm(4, size);
    int index = 0;
    for(std::map<int, MatrixXd >::iterator it = dgms.begin(); it !=
        dgms.end(); ++it){
      MatrixXd &dg = it->second;
      for(unsigned int i=0; i< dg.cols(); i++, index++){
        dgm(0, index) = dg(0, i);
        dgm(1, index) = dg(1, i);
        dgm(2, index) = it->first;
        dgm(3, index) = dg(2, i);
      }
    }
    int fSize = mrips.getMaximumFiltrationSize();

    SEXP list;
    PROTECT( list = Rf_allocVector(VECSXP, 2) );

    SEXP Rdgm;
    PROTECT( Rdgm = Rf_allocMatrix(REALSXP, dgm.rows(), dgm.cols()) );
    memcpy( REAL(Rdgm), dgm.data(), dgm.cols()*dgm.rows()*sizeof(double) );
    SET_VECTOR_ELT(list, 0, Rdgm);

    SEXP RfSize;
    PROTECT(RfSize = Rf_allocVector(INTSXP, 1));
    INTEGER(RfSize)[0] = fSize;
    SET_VECTOR_ELT(list, 1, RfSize);

    UNPROTECT(3);

    return list;

  };



  SEXP multiscale_rips(SEXP Rgmra, SEXP Rd, SEXP Rpercentile, SEXP RdType, SEXP Rsingle,
       SEXP RradiusFactor, SEXP Rrisky, SEXP RalgType, SEXP RdsType) {

    int d = *INTEGER(Rd);
    bool single = *INTEGER(Rsingle) != 0;
    float percentile = (float) *REAL(Rpercentile);
    float radiusFactor = (float) *REAL(RradiusFactor);
    int risky = *INTEGER(Rrisky);
    int algType = *INTEGER(RalgType);
    int dsType = *INTEGER(RdsType);

    SEXP Rgmra_pointer = Rf_getAttrib(Rgmra, Rf_install("gmra_ptr") );
    GMRATree<double> *gmra = static_cast<GMRATree<double> *>( R_ExternalPtrAddr(Rgmra_pointer) );

    int dType = *INTEGER(RdType);
    NodeDistance<double> *dist = getNodeDistance(dType);
    gmra->computeRadii(dist);
    gmra->computeLocalRadii(dist);
    GenericGMRANeighborhood<double> nh(gmra, dist);

    SEXP Rres;
    if( risky < 2 ){
        Rres = runMultiscaleRips(nh, d, percentile, single, radiusFactor, risky,
            algType, dsType);
    }
    else if( risky == 2){
        Rres = runMultiscaleRipsRisky1(nh, d, percentile, single, radiusFactor);
    }
    else {
        Rres = runMultiscaleRipsRisky2(nh, d, percentile, single, radiusFactor);
    }
    return Rres;
  };




  SEXP grid3d_persistence(SEXP Rd, SEXP Rf,SEXP Rmask) {

    int *dims = INTEGER(Rd);
    double *f =  REAL(Rf);
    int *mask =  INTEGER(Rmask);


    QuadFiltration filt;
    Grid3dNeighbors<double>::buildQuadFiltration(filt, dims[2], dims[1], dims[0], f, mask);

    QuadPersistentHomology<double>  ph;
    ph.run( filt  );

    QuadPersistentHomology<double>::MatrixXp dgm = ph.getDiagram();

    SEXP Rdgm;
    PROTECT( Rdgm = Rf_allocMatrix(REALSXP, dgm.rows(), dgm.cols()) );
    memcpy( REAL(Rdgm), dgm.data(), dgm.cols()*dgm.rows()*sizeof(double) );
    UNPROTECT(1);

    return Rdgm;
  };


  SEXP nn_rips( SEXP Rx, SEXP Rm, SEXP Rn, SEXP RmaxD, SEXP Rknn,
                  SEXP RincludeZeros ) {

    int knn = *INTEGER( Rknn );
    int maxD = *INTEGER( RmaxD );
    bool includeZeros = *INTEGER(RincludeZeros) != 0;
    int m = *INTEGER(Rm);
    int n = *INTEGER(Rn);
    double *x = REAL(Rx);

    FortranLinalg::DenseMatrix<double> X(m, n, x);

    FortranLinalg::DenseMatrix< int > nn( knn, X.N() );
    FortranLinalg::DenseMatrix< double > nnD( knn, X.N() );
    FortranLinalg::EuclideanMetric< double > metric;
    FortranLinalg::Distance< double >::computeKNN( X, nn, nnD, metric );

    FiltrationBuilder<double>::Neighbors N( X.N() );
    for(unsigned int i = 0; i<X.N(); i++){
      FiltrationBuilder<double>::NeighborMap NN;
      for(unsigned int j=0; j<nn.M(); j++){
        NN[ nn(j,i) ] =  nnD(j, i);
      }
      N[i] = NN;
    }

    FiltrationBuilder<double> rips;
    rips.run(N, maxD);


    Filtration &filt = rips.getFiltration();

    PersistentHomology<double>  ph;
    if(includeZeros){
      ph.run( filt, maxD );
    }
    else{
      ph.run( filt );
    }

    PersistentHomology<double>::MatrixXp dgm = ph.getDiagram();
    SEXP Rdgm;
    PROTECT( Rdgm = Rf_allocMatrix(REALSXP, dgm.rows(), dgm.cols()) );
    memcpy( REAL(Rdgm), dgm.data(), dgm.cols()*dgm.rows()*sizeof(double) );
    UNPROTECT(1);

    return Rdgm;
  };


}//end extern C
