//author: Samuel Gerber

#ifndef QUADSIMPLEX_H
#define QUADSIMPLEX_H





//typedef std::set<vertices> Simplex;


class QuadSimplex {


public:
  
  std::set<int> vertices;
  //std::list< QuadSimplex > faces;

  QuadSimplex() {
  };
  

  bool operator == (const QuadSimplex& other) const{
    return this->vertices == other.vertices;
  };
  

  bool operator < (const QuadSimplex& other) const{
    int s1 = this->vertices.size();
    int s2 = other.vertices.size();
    if(s1 !=  s2){
      return s1 < s2;
    }
    

    return this->vertices < other.vertices;
  };



  bool operator > (const QuadSimplex& other) const{

    int s1 = this->vertices.size();
    int s2 = other.vertices.size();
    if(s1 !=  s2){
      return s1 > s2;
    }


    return this->vertices > other.vertices;
  };


/*
  void addFace(QuadSimplex &q){
    faces.push_back(q);
  };
*/

  std::list<QuadSimplex> getFaces() const{

    
    std::list<QuadSimplex> faces;
    if(vertices.size() == 1){
      return faces;
    }

    QuadSimplex s;
    addRecursive(s, vertices.begin(), faces);
    

    /*
    std::set< QuadSimplex > tmp;
    addRecursive(*this, tmp);

    for(std::set<QuadSimplex>::iterator it=tmp.begin(); it !=tmp.end(); ++it){
      faces.push_back(*it);
    }
    */


    return faces;
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


  
  private:

   void addRecursive(QuadSimplex &s, std::set<int>::iterator it, std::list<QuadSimplex> &faces) const{
     if(s.vertices.size() < vertices.size()/2){
       while(it != vertices.end()){
         QuadSimplex s2 = s;
         s2.vertices.insert(*it);
         ++it;
         addRecursive(s2, it, faces);
       }
     }
     else{
       faces.push_back(s);
     }
   }

  /*
    void addRecursive(QuadSimplex s1, std::set<QuadSimplex> &add) const{
      if(s1.vertices.size() >  this->vertices.size() / 2){
        for(  std::set<int>::iterator it=s1.vertices.begin(); it != s1.vertices.end();
             ++it){
           QuadSimplex s2 = s1;
           s2.vertices.erase(*it);
           addRecursive(s2, add);
        }
      }
      else{
        add.insert(s1);
      }

    };
    */


};




#endif 

