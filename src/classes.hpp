#ifndef CLASSES_H_
#define CLASSES_H_

#include <string>

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
    std::string material;
    int material_int;
    double tmtmrcut;
    double retmrcut;
    bool tracking;
    double tt_factor;
    double rt_factor;
    double rt_exchange_constant;
    double zrconcentration;
    bool domainwall;
    bool centrepin;
    vec_t dw_dim;
    bool zrdoping;
    bool tidope;
    int ticoncentration;
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
    bool fe;                    // is it iron or not
};

// class to hold interaction information
struct int_t
{
    int iid;                // interaction id
    atom_t i, j;            // atom ids of i and j
    vec_t disp;             // unitcell displacement
    double exchange;        // exchange energy associated with atom pair
    double rij;             // atom separation
};

struct material_t
{
   std::string name;        // element name
   int id;                  // element id
};

#endif
