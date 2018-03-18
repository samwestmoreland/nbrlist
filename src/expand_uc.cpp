#include<iostream>
#include<fstream>
#include<vector>

#include "classes.hpp"
#include "data.hpp"
#include "io.hpp"

int expand_unitcell_and_substitute_zr_atoms(double zrconcentration) {

   std::cout << "expanding unitcell and substituting zr atoms.\n";

   /* reset atom count array */
   for (int i=0; i<material_specific_atom_count.size(); ++i)
      material_specific_atom_count[i] = 0;

   /* add one more array element for Zr */
   material_specific_atom_count.push_back(0);

   material_t new_material;
   new_material.name = "Zr";
   new_material.id = materials.size()-1;

   materials.push_back(new_material);
   /**/

   /* dimensions of expanded unit cell */
   exp_ucd.x = ucd.x * 3.0;
   exp_ucd.y = ucd.y * 3.0;
   exp_ucd.z = ucd.z * 3.0;

   int atom_counter = 0;

   /* replicate unit cell in 3-dimensions */
   for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
         for (int k=0; k<3; ++k) {

            for (int atom=0; atom<unitcell.size(); ++atom) {

               atom_t temp;

               temp.element = unitcell[atom].element;

               /* this is the only thing that should be different from the (nonexpanded) unitcell */
               temp.pos.x = unitcell[atom].pos.x + i * ucd.x;
               temp.pos.y = unitcell[atom].pos.y + j * ucd.y;
               temp.pos.z = unitcell[atom].pos.z + k * ucd.z;

               temp.uc.x = i;
               temp.uc.y = j;
               temp.uc.z = k;

               /* no need to preserve atom ids because this will act as the new unit cell */
               temp.aid = atom_counter;
               atom_counter ++;

               temp.hcat = unitcell[atom].hcat;

               temp.mat = unitcell[atom].mat;
               temp.fe = unitcell[atom].fe;

               /* substitute Sm with Zr if random number is below Zr concentration */
               if (temp.element == "Sm") {

                  if ((rand() % 100) < zrconcentration*100) {
                     temp.element = "Zr";
                     temp.mat = materials.size()-1;
                  }

               }

               /* add one to material specific atom count */
               material_specific_atom_count[temp.mat] ++;

               outfile
                  << temp.aid << "\t"
                  << temp.pos.x/(exp_ucd.x) << "\t"
                  << temp.pos.y/(exp_ucd.y) << "\t"
                  << temp.pos.z/(exp_ucd.z) << "\t"
                  << temp.mat << "\t"
                  << "0\t0\n";

               expanded_uc.push_back(temp);
            }
         }
      }
   }

   std::cout << "\nexpanded material array:\n\n";
   output_materials(materials, material_specific_atom_count, expanded_uc.size());

   std::cout << "expanded unitcell contains " << expanded_uc.size() << " atoms" << std::endl;


   return EXIT_SUCCESS;
}

void populate_supercell_using_expanded_unitcell() {

   /* loop through dimensions */
   for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
         for (int k=0; k<3; k++)

            /* loop through atoms in unitcell */
            for (int atom=0; atom<expanded_uc.size(); atom++) {

               atom_t temp;
               vec_t uc;
               uc.x = i;
               uc.y = j;
               uc.z = k;

               temp.aid = expanded_uc[atom].aid;
               temp.element = expanded_uc[atom].element;
               temp.mat = expanded_uc[atom].mat;

               /* replicate unitcell atoms */
               temp.pos = expanded_uc[atom].pos + (uc * exp_ucd);

               /* label unitcell coordinates */
               temp.uc = uc;

               /* dummy values for unneeded struct elements */
               temp.hcat = 0;
               temp.gid = 0;   // this isn't needed yet as the calculations are generated for a large system later

               /* place atom in array */
               expanded_sc.push_back(temp);
            }

   std::cout
      << "atoms in super cell: "
      << expanded_sc.size() << std::endl;

   array_to_rasmol(expanded_sc, "expanded_supercell");

}
