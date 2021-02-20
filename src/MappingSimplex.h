//author: Samuel Gerber

#ifndef MAPPINGSIMPLEX_H
#define MAPPINGSIMPLEX_H


// Duplicated code from Simplex.h
// Easier not to deal with Polymorphism

#include <set>
#include <list>




//typedef std::set<vertices> Simplex;


class MappingSimplex {

private:
  int scale;
  bool mapped;

public:

  std::set<int> vertices;

  MappingSimplex(int s = 0, bool m=false) :  scale(s), mapped(m){


  };



  bool operator == (const MappingSimplex& other) const{
    return this->vertices == other.vertices;
  };


  bool operator < (const MappingSimplex& other) const{
    int s1 = this->vertices.size();
    int s2 = other.vertices.size();
    if(s1 !=  s2){
      return s1 < s2;
    }

    return this->vertices < other.vertices;
  };



  bool operator > (const MappingSimplex& other) const{

    int s1 = this->vertices.size();
    int s2 = other.vertices.size();
    if(s1 !=  s2){
      return s1 > s2;
    }

    return this->vertices > other.vertices;
  };




  std::list<MappingSimplex> getFaces() const{
    std::list<MappingSimplex> faces;
    if(vertices.size() == 1){
      return faces;
    }

    for(std::set<int>::iterator it=vertices.begin(); it != vertices.end();
        ++it){
      MappingSimplex s = *this;
      s.vertices.erase(*it);
      faces.push_back(s);
    }
    return faces;
  };



  int getScale() const{
    return scale;
  };


  void setScale(int s){
    if(s > scale){
#ifdef VERBOSE
      std::cout << "increased scale" << std::endl;
#endif
    }
    scale = s;
  };

  void setMapped(bool m){
    mapped = m;
  };

  bool isMapped() const{
    return mapped;
  };


  void print() const{

#ifdef VERBOSE
    for(std::set<int>::iterator it = vertices.begin(); it !=
           vertices.end(); ++it){
      std::cout << *it << ",";
    }
    std::cout << std::endl;
#endif
  };






};




#endif

