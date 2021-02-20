//author: Samuel Gerber

#ifndef GRID3DNEIGHBORS_H
#define GRID3DNEIGHBORS_H

#include <cmath>

#include "FiltrationBuilder.h"
#include "QuadFiltration.h"

template<typename TPrecision>
class Grid3dNeighbors{
  public:
    typedef typename FiltrationBuilder<TPrecision>::NeighborMap NeighborMap;
    typedef typename NeighborMap::iterator NeighborMapIterator;
    typedef typename FiltrationBuilder<TPrecision>::Neighbors Neighbors; 
    


  public:


    Grid3dNeighbors(){ 

    };


    ~Grid3dNeighbors(){

    };



    // Build a filtarion on a grid using 8 neighbors and the difference in function values in f 
    // as the filtration entry time
    // Expects f to be linear index by z then y and then x,m i.e.
    // f(x,y,z) = f[ z + y*nz + x*(ny*nz) ]
    static void build( Neighbors &N, int nx, int ny, int nz, TPrecision *f, int *mask ){
      N.clear();
      N.reserve(nx*ny*nz);
      double maxDiff = 0;
      int slice = ny*nz;
      for(int i=0; i<nx; i++){
        int xoff = i * (slice);
        for(int j=0; j<ny; j++){
         int yoff = xoff + j*nz; 
         for(int k=0; k<nz; k++){
           int zoff = yoff + k;
           if( mask[ zoff ]  == 0 ){
             continue;
           }
           TPrecision fc = f[zoff];

           NeighborMap NN;
           // 1 - Neighbors
           for(int ii=-1; ii<=1; ii++){
             if( i+ii >= 0 && i+ii < nx){
               int xnoff = (i+ii)*(slice);
               for(int jj=-1; jj<=1; jj++){
                 if( j+jj >= 0 && j+jj < ny){
                   int ynoff = xnoff + (j+jj)*(nz);
                   for(int kk=-1; kk<=1; kk++){
                     if( k+kk >= 0 && k+kk < nz){
                       int znoff = ynoff + k + kk;
                       if(znoff != zoff){
                         TPrecision fn = f[znoff];
                         NN[znoff] = std::max(fc, fn);
                         maxDiff = std::max(maxDiff, std::fabs(fc-fn) );
                       }
                     }
                   }
                 }
               }
             }
           }
           //end computing neighbors
             
           N.push_back(NN);
         }          
        }
      } 
    };




    //Build a quad filtration (should work too?) 
    static void buildQuadFiltration(QuadFiltration &filt, int nx, int ny, int nz, TPrecision *f, int *mask ){
      
      filt.clear();
      filt.reserve( nx * ny * nz * (8 + 12 + 6 + 1) );
      int slice = ny*nz;
      for(int i=0; i<nx-1; i++){
        int x1 = i * (slice);
        int x2 = (i+1) * (slice);
        for(int j=0; j < (ny-1); j++){
          int x1y1 = x1 + j*nz; 
          int x1y2 = x1 + (j+1) * nz; 
          int x2y1 = x2 + j * nz; 
          int x2y2 = x2 + (j+1) * nz; 
          for(int k=0; k < (nz-1); k++){
            int x1y1z1 = x1y1 + k; 
            int x1y1z2 = x1y1 + k+1; 
            int x1y2z1 = x1y2 + k; 
            int x1y2z2 = x1y2 + k+1; 
            int x2y1z1 = x2y1 + k; 
            int x2y1z2 = x2y1 + k+1; 
            int x2y2z1 = x2y2 + k; 
            int x2y2z2 = x2y2 + k+1; 
           
            if( mask[ x1y1z1 ]  == 0 ){
               continue;
            }
           
            TPrecision fx1y1z1  = f[ x1y1z1 ];
            TPrecision fx1y1z2  = f[ x1y1z2 ];
            TPrecision fx1y2z1  = f[ x1y2z1 ];
            TPrecision fx1y2z2  = f[ x1y2z2 ];
            TPrecision fx2y1z1  = f[ x2y1z1 ];
            TPrecision fx2y1z2  = f[ x2y1z2 ];
            TPrecision fx2y2z1  = f[ x2y2z1 ];
            TPrecision fx2y2z2  = f[ x2y2z2 ];
           

            //Add vertices
            QuadSimplex sx1y1z1;
            sx1y1z1.vertices.insert( x1y1z1 );
            QuadSimplex sx1y1z2;
            sx1y1z2.vertices.insert( x1y1z2 );
            QuadSimplex sx1y2z1;
            sx1y2z1.vertices.insert( x1y2z1 );
            QuadSimplex sx1y2z2;
            sx1y2z2.vertices.insert( x1y2z2 );
            QuadSimplex sx2y1z1;
            sx2y1z1.vertices.insert( x2y1z1 );
            QuadSimplex sx2y1z2;
            sx2y1z2.vertices.insert( x2y1z2 );
            QuadSimplex sx2y2z1;
            sx2y2z1.vertices.insert( x2y2z1 );
            QuadSimplex sx2y2z2;
            sx2y2z2.vertices.insert( x2y2z2 );

            filt.push_back( QuadFiltrationEntry( sx1y1z1, fx1y1z1 ) );
            filt.push_back( QuadFiltrationEntry( sx1y1z2, fx1y1z2 ) );
            filt.push_back( QuadFiltrationEntry( sx1y2z1, fx1y2z1 ) );
            filt.push_back( QuadFiltrationEntry( sx1y2z2, fx1y2z2 ) );
            filt.push_back( QuadFiltrationEntry( sx2y1z1, fx2y1z1 ) );
            filt.push_back( QuadFiltrationEntry( sx2y1z2, fx2y1z2 ) );
            filt.push_back( QuadFiltrationEntry( sx2y2z1, fx2y2z1 ) );
            filt.push_back( QuadFiltrationEntry( sx2y2z2, fx2y2z2 ) );

            //Add edges
            QuadSimplex e1;
            e1.vertices.insert( x1y1z1 );
            e1.vertices.insert( x1y1z2 );
            TPrecision fe1 = std::max( fx1y1z1, fx1y1z2) ; 
            filt.push_back( QuadFiltrationEntry( e1, fe1 ) );

            QuadSimplex e2;
            e2.vertices.insert( x1y2z1 );
            e2.vertices.insert( x1y2z2 );
            TPrecision fe2 = std::max( fx1y2z1, fx1y2z2 );
            filt.push_back( QuadFiltrationEntry( e2, fe2 ) );

            QuadSimplex e3;
            e3.vertices.insert( x2y1z1 );
            e3.vertices.insert( x2y1z2 );
            TPrecision fe3 = std::max( fx2y1z1, fx2y1z2 );
            filt.push_back( QuadFiltrationEntry( e3, fe3 ) );

            QuadSimplex e4;
            e4.vertices.insert( x2y2z1 );
            e4.vertices.insert( x2y2z2 );
            TPrecision fe4 = std::max( fx2y2z1, fx2y2z2 );
            filt.push_back( QuadFiltrationEntry( e4, fe4 ) );


            QuadSimplex e5;
            e5.vertices.insert( x1y1z1 );
            e5.vertices.insert( x1y2z1 );
            TPrecision fe5 = std::max( fx1y1z1, fx1y2z1 ); 
            filt.push_back( QuadFiltrationEntry( e5, fe5 ) );

            QuadSimplex e6;
            e6.vertices.insert( x1y1z2 );
            e6.vertices.insert( x1y2z2 );
            TPrecision fe6 = std::max( fx1y1z2, fx1y2z2 );
            filt.push_back( QuadFiltrationEntry( e6, fe6 ) );

            QuadSimplex e7;
            e7.vertices.insert( x2y1z1 );
            e7.vertices.insert( x2y2z1 );
            TPrecision fe7 = std::max( fx2y1z1, fx2y2z1 );
            filt.push_back( QuadFiltrationEntry( e7, fe7 ) );

            QuadSimplex e8;
            e8.vertices.insert( x2y1z2 );
            e8.vertices.insert( x2y2z2 );
            TPrecision fe8 = std::max( fx2y1z2, fx2y2z2 );
            filt.push_back( QuadFiltrationEntry( e8, fe8 ) );

             
            QuadSimplex e9;
            e9.vertices.insert( x1y1z1 );
            e9.vertices.insert( x2y1z1 );
            TPrecision fe9 = std::max( fx1y1z1, fx2y1z1 ); 
            filt.push_back( QuadFiltrationEntry( e9, fe9 ) );

            QuadSimplex e10;
            e10.vertices.insert( x1y1z2 );
            e10.vertices.insert( x2y1z2 );
            TPrecision fe10 = std::max( fx1y1z2, fx2y1z2 );
            filt.push_back( QuadFiltrationEntry( e10, fe10 ) );

            QuadSimplex e11;
            e11.vertices.insert( x1y2z1 );
            e11.vertices.insert( x2y2z1 );
            TPrecision fe11 = std::max( fx1y2z1, fx2y2z1);
            filt.push_back( QuadFiltrationEntry( e11, fe11 ) );

            QuadSimplex e12;
            e12.vertices.insert(x1y2z2);
            e12.vertices.insert(x2y2z2);
            TPrecision fe12 = std::max( fx1y2z2, fx2y2z2);
            filt.push_back( QuadFiltrationEntry( e12, fe12 ) );


            //Add qauds
            QuadSimplex q1;
            q1.vertices.insert( x1y1z1 );
            q1.vertices.insert( x1y1z2 );
            q1.vertices.insert( x1y2z1 );
            q1.vertices.insert( x1y2z2 );
            TPrecision fq1 = std::max( std::max( fx1y1z1, fx1y1z2 ),
                                       std::max( fx1y2z1, fx1y2z2 ) );
            filt.push_back( QuadFiltrationEntry( q1, fq1 ) );

            //std::cout << fq1 << ": " << std::endl;
            //q1.print();

            QuadSimplex q2;
            q2.vertices.insert( x1y1z1 );
            q2.vertices.insert( x1y1z2 );
            q2.vertices.insert( x2y1z1 );
            q2.vertices.insert( x2y1z2 );
            TPrecision fq2 = std::max( std::max( fx1y1z1, fx1y1z2 ),
                                       std::max( fx2y1z1, fx2y1z2 ) );
            filt.push_back( QuadFiltrationEntry( q2, fq2 ) );

            QuadSimplex q3;
            q3.vertices.insert( x1y2z1 );
            q3.vertices.insert( x1y2z2 );
            q3.vertices.insert( x2y2z1 );
            q3.vertices.insert( x2y2z2 );
            TPrecision fq3 = std::max( std::max( fx1y2z1, fx1y2z2 ),
                                       std::max( fx2y2z1, fx2y2z2 ) );
            filt.push_back( QuadFiltrationEntry( q3, fq3 ) );
          

            QuadSimplex q4;
            q4.vertices.insert( x2y1z1 );
            q4.vertices.insert( x2y1z2 );
            q4.vertices.insert( x2y2z1 );
            q4.vertices.insert( x2y2z2 );
            TPrecision fq4 = std::max( std::max( fx2y1z1, fx2y1z2 ),
                                       std::max( fx2y2z1, fx2y2z2 ) );
            filt.push_back( QuadFiltrationEntry( q4, fq4 ) );
            
            //std::cout << fq4 << ": " << std::endl;
            //q4.print();

            QuadSimplex q5;
            q5.vertices.insert( x1y1z1 );
            q5.vertices.insert( x1y2z1 );
            q5.vertices.insert( x2y1z1 );
            q5.vertices.insert( x2y2z1 );
            TPrecision fq5 = std::max( std::max( fx1y1z1, fx1y2z1 ),
                                       std::max( fx2y1z1, fx2y2z1 ) );
            filt.push_back( QuadFiltrationEntry( q5, fq5 ) );

            QuadSimplex q6;
            q6.vertices.insert( x1y1z2 );
            q6.vertices.insert( x1y2z2 );
            q6.vertices.insert( x2y1z2 );
            q6.vertices.insert( x2y2z2 );
            TPrecision fq6 = std::max( std::max( fx1y1z2, fx1y2z2 ),
                                       std::max( fx2y1z2, fx2y2z2 ) );
            filt.push_back( QuadFiltrationEntry( q6, fq6 ) );
 

            //Add volume
            QuadSimplex v;
            v.vertices.insert( x1y1z1 );
            v.vertices.insert( x1y1z2 );
            v.vertices.insert( x1y2z1 );
            v.vertices.insert( x1y2z2 );
            v.vertices.insert( x2y1z1 );
            v.vertices.insert( x2y1z2 );
            v.vertices.insert( x2y2z1 );
            v.vertices.insert( x2y2z2 );
            filt.push_back( QuadFiltrationEntry( v, std::max(fq5, fq6) ) );

          }          
        }
      }

      std::cout << "QuadFiltration size: " << filt.size() << std::endl;
      std::sort( filt.begin(), filt.end() );
      QuadFiltrationIterator it = std::unique( filt.begin(), filt.end() );
      filt.erase( it, filt.end() );
      std::cout << "QuadFiltration unique size: " << filt.size() << std::endl;

    };

};

#endif 

