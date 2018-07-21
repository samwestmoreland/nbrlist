#ifndef CLASSES_H_
#define CLASSES_H_

#include <string>
#include <cmath>

class lattparameters_t {
   public:
      int concentration;
      double a;
      double b;
      double c;
};

class vec_t {
public:
    double x;
    double y;
    double z;

    // overload '-' operator to subtract two vec_t objects
    vec_t operator-(const vec_t& v)
    {
        vec_t vec;
        vec.x = this->x - v.x;
        vec.y = this->y - v.y;
        vec.z = this->z - v.z;
        return vec;
    }

    vec_t operator*(const vec_t& v)
    {
        vec_t vec;
        vec.x = this->x * v.x;
        vec.y = this->y * v.y;
        vec.z = this->z * v.z;
        return vec;
    }

    vec_t operator+(const vec_t& v)
    {
        vec_t vec;
        vec.x = this->x + v.x;
        vec.y = this->y + v.y;
        vec.z = this->z + v.z;
        return vec;
    }
};

struct parameter_t {

    double rcut_tt;
    double rcut_rt;
    int rt_shell;
    double tt_factor;
    double rt_factor;
    double rt_constant;

    double zr_concentration;
    int ti;    // ti concentration

    bool tracking;
    bool domainwall;
    bool centrepin;

    bool rt_shell_cut;

    vec_t dw_dim;
};

// class for all atoms in system
class atom_t {
   public:
      std::string element;        // element string
      // vec_t spin;                 // spin components
      vec_t pos;                  // coordinates
      vec_t uc;                   // unitcell coordinates
      int aid;                    // atom id (unique within unitcell)
      int gid;                    // global atom id
      int hcat;                   // height category
      int mat;                    // material

      bool is_tm() {

         if (element == "Fe" ||
             element == "Fe8i" ||
             element == "Fe8j" ||
             element == "Fe8f" ||
             element == "Co" ||
             element == "Co8i" ||
             element == "Co8j" ||
             element == "Co8f" ||
             element == "Fe16k1" ||
             element == "Fe16k2" ||
             element == "Fe8j1" ||
             element == "Fe8j2" ||
             element == "Fe4e" ||
             element == "Fe4c")
            return true;
         else return false;
      }

      bool is_re() {

         if (element == "Nd" ||
             element == "Sm" ||
             element == "Nd4f" ||
             element == "Nd4g")
            return true;
         else return false;
      }
};

class material_t {
   public:
      std::string name;        // element name
      vec_t ucd;               // unitcell dimensions

      int id() {
         if (name == "bccfe") return 1;
         else if (name == "ndfeb") return 2;
         else if (name == "ndfe12") return 3;
         else if (name == "smfe12") return 4;
         else if (name == "smco12") return 5;
         else if (name == "interface") return 6;
         else if (name == "ndfeti12") return 7;
         else return 0;
      }
};

class pair_t {
   public:
      int iid;             // interaction id
      atom_t i, j;         // atoms i and j
      vec_t ucd;           // unitcell displacement
      double exchange;     // exchange energy associated with interaction

      /* check atoms aren't interacting with themselves */
      bool are_same_atom() {
         if (i.gid == j.gid)
            return true;
         else return false;
      }

      /* calculate inter-atomic separation */
      double rij() {
         double dx = j.pos.x - i.pos.x;
         double dy = j.pos.y - i.pos.y;
         double dz = j.pos.z - i.pos.z;

         return sqrt(dx*dx+dy*dy+dz*dz);
      }

      /* are the atoms within the cut off radii? */
      bool are_within_range(double rcut_rt, double rcut_tt) {
         if ((i.is_re() && j.is_tm()) || (i.is_tm() && j.is_re())) {
            if (rij() <= rcut_rt) return true;
            else return false;
         }

         else if (i.is_tm() && j.is_tm()) {
            if (rij() <= rcut_tt) return true;
            else return false;
         }

         else return false;
      }
};

#endif
