// scaling of m for generic bcc system

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>

class vec_t
{
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

// class for all atoms in system
class atom_t
{
public:
    std::string element;        // element string
    vec_t spin;                 // spin components
    vec_t pos;                  // coordinates
    vec_t uc;                   // unitcell coordinates
    int aid;                    // atom id (unique within unitcell)
    int gid;                    // global atom id
    int hcat;                   // height category
    int mat;                    // material
};

// class to hold interaction information
class int_t
{
public:
    int iid;                // interaction id
    atom_t i, j;            // atom ids of i and j
    vec_t disp;             // unitcell displacement
    double exchange;        // exchange energy associated with atom pair
};

double calc_rij (vec_t i, vec_t j);     // distance calculation. pass positions

// exchange function
double calc_jij (std::string i_type, std::string j_type, double rij) {

    double A=36.9434;
    double B=1.25094;
    double C=-0.229572;
    const double Fe_ratio_ndfeb=0.69618016759*1.07692307692; // 560/520 = 1.07692307692
    const double Fe_ratio_ndfe12=1.15*1.07692307692; // 560/520 = 1.07692307692
    // ^^ ndfeb

    const double J0Nd_ndfeb=Fe_ratio_ndfeb*4.06835e-20/16.0;
    const double J0Nd_ndfe12=Fe_ratio_ndfeb*4.06835e-20/16.0;

    // Nd-Nd
    if(i_type=="Nd" && j_type=="Nd"){
        return 0.0;
    }

    // Nd-Fe (step function at r = 4A)
    if((i_type=="Fe" && j_type=="Nd") || (i_type=="Nd" && j_type=="Fe"))
    {
        if(rij<=4.0) return 0.33*J0Nd_ndfeb;       // ndfeb
        else return 0.0;
    }

    // Fe-Fe (cutoff at r = 5.74A)
    else if(i_type=="Fe" && j_type=="Fe"){
        //if(rij<=5.0) return -2.0*2.179872e-21*(A*exp(-B*rij)+C);
        // Correct for Tc = 600
        if(rij<=5.0) return 2.0*2.179872e-21*(A*exp(-B*rij)+C)*Fe_ratio_ndfeb;
        else return 0.0;

        //        if (rij <= 2.6) return 3.117407e-21;
        //        else if (rij <= 3.0) return 1.7742e-21;
        //        else return 0.0;
    }

    // B-x
    else return 0.0;

}

// function to calculate distance
double calc_rij (vec_t i, vec_t j)
{
    vec_t d = j-i;
    double distance = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);

    return distance;
}

int main (int argc, char *argv[]) {

    if (argc!=3)
    {
        std::cout << "insert arguments: input file (unscaled), tracked cell dimension" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    double ntr = atof(argv[2]); /* number of tracked cells */

    // =============================
    // store unitcell atoms in array
    // =============================

    // open unitcell coordinates file      // change
    std::ifstream uccoords (argv[1]);

    // initialise array to hold unitcell info
    std::vector<atom_t> unitcell;

    /*
    *  lattice parameters  // change
    */

    // double ucd[3] = {8.701,8.701,4.844}; // --- ndfe12

    // ******** ndfeb
    // double a = 8.8;
    // double c = 12.2;

    // ******** ndfe12
    // double a = 8.701;
    // double c = 4.844;

    double a = 1;
    double c = 1;

    vec_t ucd;
    ucd.x = a;
    ucd.y = a;
    ucd.z = c;

    std::vector<int> material_count(6, 0);

    /* read coordinates into array */
    if (uccoords)
    {
        atom_t tmp; // temporary variable
        int aid_counter = 0;

        // read file into array
        while (uccoords >> tmp.element >> tmp.pos.x >> tmp.pos.y >> tmp.pos.z) {

            // if atom magnetic, place in array
            if (!(tmp.element=="B")) {
                tmp.aid = aid_counter;

                // scale coordinates
                tmp.pos = tmp.pos * ucd;

                /* bottom Nd2Fe14B layer */
                if (tmp.pos.z < 0.0)
                {
                    if (tmp.element == "Nd")
                    {
                        if (!(tmp.pos.z < -0.0 && tmp.pos.z > -6.0))
                        {
                            tmp.mat = 0;
                            ++material_count.at(0);
                        }
                        else
                        {
                            tmp.mat = 5; /* first Nd plane */
                            ++material_count.at(5);
                        }
                    }
                    else if (tmp.element == "Fe")
                    {
                        tmp.mat = 1;
                        ++material_count.at(1);
                    }
                }

                /* bcc Fe layer */
                if (tmp.pos.z > 0.0 && tmp.pos.z <= 109.0)
                {
                    tmp.mat = 2;
                    ++material_count.at(2);
                }

                /* top Nd2Fe14B layer */
                if (tmp.pos.z > 109.0)
                {
                    if (tmp.element == "Nd")
                    {
                        tmp.mat = 3;
                        ++material_count.at(3);
                    }

                    else if (tmp.element == "Fe")
                    {
                        tmp.mat = 4;
                        ++material_count.at(4);
                    }
                }
unitcell.push_back(tmp);
                ++aid_counter;
            }
        }

        std::cout << "coordinates read.\n\n";
    }

    else std::cerr << "no input file found\n";

    for (int i = 0; i < material_count.size(); ++i)
        std::cout << "atoms of material " << i << ": "
                  << material_count.at(i) << std::endl;

    std::cout << "\nnumber of atoms in unitcell: " << unitcell.size() << std::endl;

    /*
     * replicate unitcell to make supercell
     */

    std::vector<atom_t> super;
    super.reserve(unitcell.size()*3*3*3);

    // loop through dimensions
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    for (int k=0; k<1; k++)

    /* loop through atoms in unitcell */
    for (int atom=0; atom<unitcell.size(); atom++)
    {
        atom_t tmp;
        vec_t uc;
        uc.x = i;
        uc.y = j;
        uc.z = k;

        tmp.aid = unitcell[atom].aid;
        tmp.element = unitcell[atom].element;
        tmp.mat = unitcell[atom].mat;

        /* replicate unitcell atoms */
        tmp.pos = unitcell[atom].pos + (uc * ucd);

        /* label unitcell coordinates */
        tmp.uc = uc;

        /* place atom in array */
        super.push_back(tmp);
    }

    /* open file for rasmol output */
    std::ofstream sxyz ("super.xyz");

    // if file open
    if (sxyz)
    {
        sxyz << super.size() << "\n\n";

        for (int i=0; i<super.size(); i++)
            sxyz << super[i].element << "\t"
            << super[i].pos.x   << "\t"
            << super[i].pos.y   << "\t"
            << super[i].pos.z   << "\t"
            << super[i].uc.x    << "\t"
            << super[i].uc.y    << "\t"
            << super[i].uc.z    << "\n";

        sxyz.close();
        std::cout << "super.xyz file generated.\n";
    }

    else std::cerr << "*** could not open super.xyz for output ***" << std::endl;
    std::cout << "number of atoms in supercell: " << super.size() << std::endl;

    /*
    * neighbour list
    */

    std::vector<int_t> ucints;

    // cut off radius in angstroms
    double rcut = 3.0;

    // central cell location
    int start = (super.size()-unitcell.size())/2;
    int end   = (super.size()+unitcell.size())/2;

    int count = 0;

    // loop through central cell
    for (int i=start; i<end; ++i)
    {
        // loop through super cell (looking for atom j)
        for (int j=0; j<super.size(); ++j)
        {
            /* calculate interatomic distance */
            double rij = calc_rij(super[i].pos, super[j].pos);

            /* if distance less than rcut and not same atom */
            if (rij<rcut && rij!=0)
            {
                // create a pair
                int_t tmp;
                tmp.i = super[i];
                tmp.j = super[j];

                // calculate unitcell displacement
                tmp.disp = tmp.j.uc - tmp.i.uc;

                // calculate exchange energy
                tmp.exchange = calc_jij(tmp.i.element, tmp.j.element, rij);
                if (tmp.exchange!=0) ucints.push_back(tmp);
            }
        }
    }

    //        singlecell << interactions[i].iid << "\t"
    //            << interactions[i].i.gid << "\t"
    //            << interactions[i].j.gid << "\t"
    //            << interactions[i].disp.x << "\t"
    //            << interactions[i].disp.y << "\t"
    //            << interactions[i].disp.z << "\t"
    //            << interactions[i].exchange << "\n";

    std::cout << "number of interactions in unit cell: " << ucints.size() << std::endl;

    /*
    * generate large system
    */

    // system dimensions
    vec_t sd;
    sd.x = 1;
    sd.y = 1;
    sd.z = 1;

    std::cout << "\ndimensions of large system [A]: ("
    << sd.x * ucd.x << ", "
    << sd.y * ucd.y << ", "
    << sd.z * ucd.z << ")\n";

    int ns = sd.x * sd.y * sd.z * unitcell.size(); // number of atoms in system
    std::cout << "number of atoms in system: " << ns << "\n";

    // open file for rasmol output of system
    std::ofstream sysxyz ("system.xyz");
    sysxyz << ns << "\n\n";

    // open file for ucf output
    std::ofstream ucf ("ndfeb.ucf");            // change

    ucf << "# Unit cell size:\n"
    <<  sd.x*ucd.x << "\t"
    <<  sd.y*ucd.y << "\t"
    <<  sd.z*ucd.z << "\n"
    << "# Unit cell vectors:\n"
    << "1.0  0.0  0.0\n"
    << "0.0  1.0  0.0\n"
    << "0.0  0.0  1.0\n"
    << "# Atoms num, id cx cy cz mat lc hc\n"
    << ns << "\n";

    // array to hold all atoms in system
    std::vector < std::vector < std::vector < std::vector <atom_t> > > > system;

    // global id counter
    int gid_counter = 0;
    int tracked = 0;

    /* resize vectors */
    system.resize(sd.x);
    for (int i=0; i<sd.x; i++) {
        system[i].resize(sd.y);
        for (int j=0; j<sd.y; j++) {
            system[i][j].resize(sd.z);
            for (int k=0; k<sd.z; k++) {

                /* loop through unitcell atoms */
                for (int atom=0; atom<unitcell.size(); ++atom) {

                    atom_t tmp;
                    vec_t uc;
                    uc.x = i;
                    uc.y = j;
                    uc.z = k;

                    tmp.aid = unitcell[atom].aid;
                    tmp.gid = gid_counter;

                    tmp.element = unitcell[atom].element;

                    tmp.mat = unitcell[atom].mat;

                    /*
                    * cell tracking
                    */

                    // if ( (i<((sd.x+ntr/2)-sd.x/2))
                    //    &&(i>((sd.x-ntr/2)-sd.x/2))
                    //    &&(j<((sd.y+ntr/2)-sd.y/2))
                    //    &&(j>((sd.y-ntr/2)-sd.y/2))
                    //    &&(k<((sd.z+ntr/2)-sd.z/2))
                    //    &&(k>((sd.z-ntr/2)-sd.z/2)))
                    // {
                    //     tmp.mat = 1;
                    //     tmp.element = 'h';
                    //     ++tracked;
                    // }

                    tmp.pos = unitcell[atom].pos + uc*ucd;
                    ++gid_counter;

                    sysxyz << tmp.element << "\t"
                    << tmp.pos.x << "\t"
                    << tmp.pos.y << "\t"
                    << tmp.pos.z << "\n";

                    // normalise coordinates to large system dimensions
                    double sysx = tmp.pos.x / double(ucd.x) / double(sd.x);
                    double sysy = tmp.pos.y / double(ucd.y) / double(sd.y);
                    double sysz = tmp.pos.z / double(ucd.z) / double(sd.z);

                    // height categorisation
                    double hinc = 6;            // change this
                    tmp.hcat = sysz*hinc;

                    // output to file
                    ucf << tmp.gid << "\t"
                    << sysx << "\t"
                    << sysy << "\t"
                    << sysz << "\t"
                    << tmp.mat << "\t"
                    << 0 << "\t"
                    << tmp.hcat << "\n";

                    system[i][j][k].push_back(tmp);

                }
            }
        }
    }

    std::cout << "number of tracked cells: " << tracked/unitcell.size() << "\n";

    sysxyz.close();

    ucf << "# Interactions n exctype, id i j dx dy dz Jij\n";

    // ===============================================
    // determine interactions for every atom in system
    // using pre-calculated neighbour list
    // ===============================================

    // initialise array to hold interactions for whole system
    std::vector<int_t> interactions;

    // interaction counter
    int int_counter = 0;

    // loop through every atom in system
    for (int i=0; i<sd.x; i++)
    for (int j=0; j<sd.y; j++)
    for (int k=0; k<sd.z; k++)

    // first atom
    for (int atom=0; atom<unitcell.size(); atom++)

    // loop through interaction information
    for (int p=0; p<ucints.size(); p++) {

        // if interaction info refers to correct atom
        if (system[i][j][k][atom].aid == ucints[p].i.aid)
        {
            int_t tmp;

            tmp.iid = int_counter;
            tmp.i.gid = system[i][j][k][atom].gid;

            tmp.i.mat = system[i][j][k][atom].mat;
            tmp.i.element = system[i][j][k][atom].element;

            tmp.j.element = ucints[p].j.element;

            tmp.exchange = ucints[p].exchange;

            // assume atom j is within system
            // to begin with
            tmp.disp.x = 0;
            tmp.disp.y = 0;
            tmp.disp.z = 0;

            // check if atom j is within system boundaries
            int ucx = i + ucints[p].disp.x;
            int ucy = j + ucints[p].disp.y;
            int ucz = k + ucints[p].disp.z;

            // only one of these conditions
            // needs to be satisfied for the
            // atom to be out of bounds
            if ( ucx < 0 || ucx >= sd.x ||
                 ucy < 0 || ucy >= sd.y ||
                 ucz < 0 || ucz >= sd.z )
                {

                    // periodic boundaries conditions
                    if (ucx<0) {
                        ucx += sd.x;
                        tmp.disp.x = -1;

                        //    // determine mat of atom j
                        //    if (tmp.i.mat==5)
                        //    {
                        //        if (tmp.j.element=="Nd") tmp.j.mat = 5;
                        //        else tmp.j.mat = 6;
                        //    }
                        //
                        //    else if (tmp.i.mat==6)
                        //    {
                        //        if (tmp.j.element=="Fe") tmp.j.mat = 6;
                        //        else tmp.j.mat = 5;
                        //    }

                    }

                    if (ucy<0) {
                        ucy += sd.y;
                        tmp.disp.y = -1;
                    }

                    if (ucz<0) {
                        ucz += sd.z;
                        tmp.disp.z = -1;
                    }
                    if (ucx>=sd.x) {
                        ucx -= sd.x;
                        tmp.disp.x = 1;
                    }
                    if (ucy>=sd.y) {
                        ucy -= sd.y;
                        tmp.disp.y = 1;
                    }
                    if (ucz>=sd.z) {
                        ucz -= sd.z;
                        tmp.disp.z = 1;
                    }

                    // having changed uc coordinates
                    // obtain j.gid
                    tmp.j.gid = system[ucx][ucy][ucz][ucints[p].j.aid].gid;

                }

                // else it is within bounds
                // so simply extract j.gid
                else tmp.j.gid = system[ucx][ucy][ucz][ucints[p].j.aid].gid;

                /*
                // if at top afm boundary
                // switch sign of exchange
                // and amplify exchange to pin wall
                if (k==sd[2])
                {
                if (tmp.dz == -1) tmp.exchange *= -100;
            }
            else if (k==sd[2]-1)
            {
            if (tmp.dz == 1) tmp.exchange *= -100;
        }
        */

        interactions.push_back(tmp);

        // increment interaction id
        ++int_counter;

    }
}

// std::cout << "number of interactions in system: " << interactions.size() << "\n\n";

ucf << interactions.size() << "\t0\n";

// output interaction info to file
for (int i=0; i<interactions.size(); i++)

ucf << interactions[i].iid << "\t"
<< interactions[i].i.gid << "\t"
<< interactions[i].j.gid << "\t"
<< interactions[i].disp.x << "\t"
<< interactions[i].disp.y << "\t"
<< interactions[i].disp.z << "\t"
<< interactions[i].exchange << "\n";

ucf.close();

return 0;

}
