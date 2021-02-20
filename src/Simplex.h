//author: Samuel Gerber

#ifndef SIMPLEX_H
#define SIMPLEX_H


#include <set>




//typedef std::set<vertices> Simplex;


class Simplex{


public:

  std::set<int> vertices;

  Simplex() {
  };


  bool operator == (const Simplex& other) const{
    return this->vertices == other.vertices;
  };


  bool operator < (const Simplex& other) const{
    int s1 = this->vertices.size();
    int s2 = other.vertices.size();
    if(s1 !=  s2){
      return s1 < s2;
    }


    return this->vertices < other.vertices;
  };



  bool operator > (const Simplex& other) const{

    int s1 = this->vertices.size();
    int s2 = other.vertices.size();
    if(s1 !=  s2){
      return s1 > s2;
    }


    return this->vertices > other.vertices;
  };




  std::list<Simplex> getFaces() const{
    std::list<Simplex> faces;
    if(vertices.size() == 1){
      return faces;
    }

    for(std::set<int>::iterator it=vertices.begin(); it != vertices.end();
        ++it){
      Simplex s = *this;
      s.vertices.erase(*it);
      faces.push_back(s);
    }
    return faces;
  };




  void print() const{

//#ifdef VERBOSE
    for(std::set<int>::iterator it = vertices.begin(); it !=
           vertices.end(); ++it){
      std::cout << *it << ",";
    }
    std::cout << std::endl;
//#endif
  };


};
#endif

