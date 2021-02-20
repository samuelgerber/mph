//author: Samuel Gerber


#ifndef GMRAIDNODE_H
#define GMRAIDNODE_H


#include "GMRATree.h"


template <typename TPrecision>
class GMRAIDNode : public GMRANodeDecorator<TPrecision>{

  private:
    int id;

  public:
    GMRAIDNode(GMRANode<TPrecision> *n, int nodeID = -1) : GMRANodeDecorator<TPrecision>(n), id(nodeID) {
    };

    void setID(int nodeID){
      id=nodeID;
    };

    int getID(){
      return id;
    }
};



template <typename TPrecision>
class GMRAIDDecorator : public Decorator<TPrecision>{

  private:
    int idCounter;
  public:
    GMRAIDDecorator() : idCounter(0) {
    };

    void resetCounter(){
      idCounter=0;
    };

    virtual GMRANode<TPrecision> *decorate(GMRANode<TPrecision> *node, GMRANode<TPrecision> *parent){

      GMRANode<TPrecision> *dNode = new GMRAIDNode<TPrecision>(node, idCounter);
      ++idCounter;
      return dNode;
    };

    int getMaxID(){
      return idCounter;
    };

};



#endif
