#include <iostream>
#include <fstream>

/* project header files */
#include "./initialise.hpp"
#include "./data.hpp"
#include "./main.hpp"
#include "./io.hpp"

/* function prototypes */
void populate_supercell();

int calculate_interactions() {

   std::cout << "\ncreating and populating super cell...\n";

   populate_supercell();

   std::cout << "calculating interactions...\n\n";

   if (sys.material_int == 2 || sys.material_int == 8)
      std::cout <<
         "TM-TM exchange factor = " << sys.tt_factor << "\n" <<
         "RE-TM exchange factor = " << sys.rt_factor << "\n";

   /* first and last index of central unitcell */
   int start;
   int end;

   int interaction_count = 0;
   int tmtm_interaction_count = 0;
   int retm_interaction_count = 0;

   double total_neighbour_distance = 0;
   double total_tmtm_neighbour_distance = 0;
   double total_retm_neighbour_distance = 0;

      start = (supercell.size()-unitcell.size())/2;
      end   = (supercell.size()+unitcell.size())/2;

      /* loop through central cell */
      for (int i=start; i<end; ++i) {

         /* loop through super cell (looking for atom j) */
         for (int j=0; j<supercell.size(); ++j) {

            /* calculate interatomic distance */
            double rij = calculate_rij(supercell[i].pos, supercell[j].pos);

            /* if distance less than rcut and not same atom */
            if (rij < sys.tmtmrcut && rij > 1e-10) {

               /* add neighbour distance to total */
               total_neighbour_distance += rij;

               /* create an interaction */
               int_t temp;
               temp.i = supercell[i];
               temp.j = supercell[j];

               temp.iid = interaction_count;

               /* unitcell displacement */
               temp.disp = temp.j.uc - temp.i.uc;

               /* calculate exchange energy */
               temp.exchange = calculate_jij(temp.i.element, temp.j.element, rij, sys.material_int);

               temp.rij = rij;

               /* put interaction into array */
               if (temp.exchange!=0) {

                  uc_interactions.push_back(temp);
                  interaction_count ++;

                  /* if interaction is fe-fe or fe-co */
                  if (
                        (temp.i.element == "Fe8i" || temp.i.element == "Fe8j" || temp.i.element == "Fe8f" || temp.i.element == "Fe" || temp.i.element == "Co") &&
                        (temp.j.element == "Fe8i" || temp.j.element == "Fe8j" || temp.j.element == "Fe8f" || temp.j.element == "Fe" || temp.j.element == "Co")
                     ) {

                     tmtm_interaction_count ++;
                     total_tmtm_neighbour_distance += rij;
                  }

                  /* if interaction is fe-re */
                  else if (
                        ((temp.i.element == "Fe8i" || temp.i.element == "Fe8j" || temp.i.element == "Fe8f" || temp.i.element == "Fe" || temp.i.element == "Co") &&
                         (temp.j.element == "Sm" || temp.j.element == "Nd"))
                        ||
                        ((temp.j.element == "Fe8i" || temp.j.element == "Fe8j" || temp.j.element == "Fe8f" || temp.j.element == "Fe" || temp.j.element == "Co") &&
                         (temp.i.element == "Sm" || temp.i.element == "Nd"))
                        ) {
                     retm_interaction_count ++;
                     total_retm_neighbour_distance += rij;
                  }
               }

            }
         }
      }

      std::cout
         << "total interactions: "
         << uc_interactions.size() << std::endl;

      std::cout
         << "tm-tm interactions: "
         << tmtm_interaction_count << std::endl;

      std::cout
         << "re-tm interactions: "
         << retm_interaction_count << std::endl;

      std::cout << std::endl;

      std::cout
         << "mean neighbour distance: "
         << total_neighbour_distance/interaction_count
         << " A" << std::endl;

      std::cout
         << "mean tm-tm neighbour distance: "
         << total_tmtm_neighbour_distance/tmtm_interaction_count
         << " A" << std::endl;

      std::cout
         << "mean re-fe neighbour distance: "
         << total_retm_neighbour_distance/retm_interaction_count
         << " A" << std::endl;

      std::cout << std::endl;

      outfile
         << "# Interactions n exctype, id i j dx dy dz Jij\n"
         << uc_interactions.size() << "\tisotropic\n";

      for (int i=0; i<uc_interactions.size(); ++i)
         outfile
            << uc_interactions[i].iid       << "\t"
            << uc_interactions[i].i.aid     << "\t"
            << uc_interactions[i].j.aid     << "\t"
            << uc_interactions[i].disp.x    << "\t"
            << uc_interactions[i].disp.y    << "\t"
            << uc_interactions[i].disp.z    << "\t"
            << uc_interactions[i].exchange  << "\n";

      std::cout << "interaction data output to \'output.ucf\'\n";

      return EXIT_SUCCESS;
}

/* calculate exchange energy */
double calculate_jij(std::string const& i_type,
                     std::string const& j_type,
                     double rij,
                     int exchange_fn) {

   switch (exchange_fn) {

      case 1 : { // bccfe
         /* parameters from Pajda (2001) */
         /* factor 2 to correct for Hamiltonian */
         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;
         return c*(a/(rij*rij*rij)-b);
      }
         break;

      case 2 : { // ndfeb

         // return jij_ndfeb(i_type, j_type, rij);

         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;

         double ndfeb_tt_factor = sys.tt_factor;
         double ndfeb_rt_factor = sys.rt_factor;

         /* Nd-Nd */
         if (i_type=="Nd" && j_type=="Nd") return 0.0;

         /* Sm-Fe */
         // if ((i_type=="Fe8i" && j_type=="Sm") ||
         //     (i_type=="Sm" && j_type=="Fe8i") ||
         //     (i_type=="Fe8j" && j_type=="Sm") ||
         //     (i_type=="Sm" && j_type=="Fe8j") ||
         //     (i_type=="Fe8f" && j_type=="Sm") ||
         //     (i_type=="Sm" && j_type=="Fe8f") ) {
         //
         //
         //    if (rij<=4.0) return smfe12_rt_factor * c*(a/(rij*rij*rij)-b);
         //    else return 0.0;
         // }


         else if ((i_type == "Nd" && j_type == "Fe") ||
                  (i_type == "Fe" && j_type == "Nd"))
         {
            return ndfeb_rt_factor * c*(a/(rij*rij*rij)-b);
         }

         // Fe-Fe (cutoff at r = 5.74A)
         else if (i_type=="Fe" && j_type=="Fe")
            return ndfeb_tt_factor * c*(a/(rij*rij*rij)-b);

         else return 0.0;

               }
               break;

      case 3 : { // ndfe12

         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;

         double ndfe12_tt_factor = sys.tt_factor;
         double ndfe12_rt_factor = sys.rt_factor;

         /* Nd-Nd */
         if(i_type=="Nd" && j_type=="Nd") return 0.0;

         /* Nd-Fe */
         if ((i_type=="Fe8i" && j_type=="Nd") ||
             (i_type=="Nd" && j_type=="Fe8i") ||
             (i_type=="Fe8j" && j_type=="Nd") ||
             (i_type=="Nd" && j_type=="Fe8j") ||
             (i_type=="Fe8f" && j_type=="Nd") ||
             (i_type=="Nd" && j_type=="Fe8f") ) {

            if (rij <= sys.retmrcut) return sys.rt_exchange_constant;
            else return 0.0;
         }

         // Fe-Fe
         else if ((i_type=="Fe8i" && j_type=="Fe8i") ||
                  (i_type=="Fe8i" && j_type=="Fe8j") ||
                  (i_type=="Fe8i" && j_type=="Fe8f") ||
                  (i_type=="Fe8j" && j_type=="Fe8i") ||
                  (i_type=="Fe8j" && j_type=="Fe8j") ||
                  (i_type=="Fe8j" && j_type=="Fe8f") ||
                  (i_type=="Fe8f" && j_type=="Fe8i") ||
                  (i_type=="Fe8f" && j_type=="Fe8j") ||
                  (i_type=="Fe8f" && j_type=="Fe8f"))

            return ndfe12_tt_factor * c*(a/(rij*rij*rij)-b);

         else return 0.0;
    	}
      	break;

      case 4 : { /* smfe12 */

         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;

         double smfe12_tt_factor = sys.tt_factor;

         /* these factors come from simulations with re-tm cut-off 4.0 A */
         // double smfe12_tt_factor = 0.65;
         // double smfe12_rt_factor = 0.20;

         /* Sm-Sm */
         if(i_type=="Sm" && j_type=="Sm") return 0.0;

         /* Sm-Fe */
         if ((i_type=="Fe8i" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8i") ||
             (i_type=="Fe8j" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8j") ||
             (i_type=="Fe8f" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8f") ) {

            if (rij <= sys.retmrcut) return sys.rt_exchange_constant;
            else return 0.0;
         }

         // Fe-Fe (cutoff at r = 5.74A)
         else if ((i_type=="Fe8i" && j_type=="Fe8i") ||
                  (i_type=="Fe8i" && j_type=="Fe8j") ||
                  (i_type=="Fe8i" && j_type=="Fe8f") ||
                  (i_type=="Fe8j" && j_type=="Fe8i") ||
                  (i_type=="Fe8j" && j_type=="Fe8j") ||
                  (i_type=="Fe8j" && j_type=="Fe8f") ||
                  (i_type=="Fe8f" && j_type=="Fe8i") ||
                  (i_type=="Fe8f" && j_type=="Fe8j") ||
                  (i_type=="Fe8f" && j_type=="Fe8f"))

            return smfe12_tt_factor * c*(a/(rij*rij*rij)-b);

         else return 0.0;
      }
         break;

      case 5 : { /* smzrfe12 */

         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;

         /* Sm-Sm */
         if(i_type=="Sm" && j_type=="Sm") return 0.0;

         /* Sm-Fe */
         if ((i_type=="Fe8i" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8i") ||
             (i_type=="Fe8j" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8j") ||
             (i_type=="Fe8f" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe8f") ) {

            if (rij<=4.0) return sys.rt_factor * c*(a/(rij*rij*rij)-b);
            else return 0.0;
         }

         // Fe-Fe (cutoff at r = 5.74A)
         else if ((i_type=="Fe8i" && j_type=="Fe8i") ||
                  (i_type=="Fe8i" && j_type=="Fe8j") ||
                  (i_type=="Fe8i" && j_type=="Fe8f") ||
                  (i_type=="Fe8j" && j_type=="Fe8i") ||
                  (i_type=="Fe8j" && j_type=="Fe8j") ||
                  (i_type=="Fe8j" && j_type=="Fe8f") ||
                  (i_type=="Fe8f" && j_type=="Fe8i") ||
                  (i_type=="Fe8f" && j_type=="Fe8j") ||
                  (i_type=="Fe8f" && j_type=="Fe8f"))

            return sys.tt_factor * c*(a/(rij*rij*rij)-b);

         else return 0.0;
      }
         break;

      case 6 : /* interface */

         return jij_ndfeb(i_type, j_type, rij);
         break;

      case 7 : /* interface mirror */

         return jij_ndfeb(i_type, j_type, rij);
         break;

      case 8 : { /* smco12 */

         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;

         double smco12_tt_factor = sys.tt_factor;
         // double smco12_rt_factor = sys.rt_factor;

         /* Sm-Sm */
         if(i_type=="Sm" && j_type=="Sm") return 0.0;

         /* Sm-Co */
         if ((i_type=="Co" && j_type=="Sm") || (i_type=="Sm" && j_type=="Co")) {
            if (rij <= sys.retmrcut) return sys.rt_exchange_constant;
            else return 0.0;
         }

         // Co-Co
         else if (i_type=="Co" && j_type=="Co")

            return smco12_tt_factor * c*(a/(rij*rij*rij)-b);

         else return 0.0;

         }
         break;

      case 9 : { /* smfeti12 */

         double a = 121.00658;
         double b = 1.72543313196278;
         double c = 1e-21;

         double smfe12_tt_factor = sys.tt_factor;
         double smfe12_rt_factor = sys.rt_factor;

         /* Sm-Sm */
         if(i_type=="Sm" && j_type=="Sm") return 0.0;

         /* Sm-Fe */
         if ((i_type=="Fe" && j_type=="Sm") ||
             (i_type=="Sm" && j_type=="Fe")) {

            return smfe12_rt_factor * c*(a/(rij*rij*rij)-b);
         }

         // Fe-Fe (cutoff at r = 5.74A)
//         else if ((i_type=="Fe8i" && j_type=="Fe8i") ||
//                  (i_type=="Fe8i" && j_type=="Fe8j") ||
//                  (i_type=="Fe8i" && j_type=="Fe8f") ||
//                  (i_type=="Fe8j" && j_type=="Fe8i") ||
//                  (i_type=="Fe8j" && j_type=="Fe8j") ||
//                  (i_type=="Fe8j" && j_type=="Fe8f") ||
//                  (i_type=="Fe8f" && j_type=="Fe8i") ||
//                  (i_type=="Fe8f" && j_type=="Fe8j") ||
//                  (i_type=="Fe8f" && j_type=="Fe8f"))

         else if ((i_type == "Fe" && j_type == "Fe"))
            return smfe12_tt_factor * c*(a/(rij*rij*rij)-b);

         else if ((i_type == "Ti" && j_type == "Ti"))
            return 0.0;

         else if ((i_type == "Ti" && j_type == "Sm"))
            return 0.0;

         else if ((i_type == "Sm" && j_type == "Ti"))
            return 0.0;

         else if ((i_type == "Ti" && j_type == "Fe"))
            return 0.0;

         else if ((i_type == "Fe" && j_type == "Ti"))
            return 0.0;

         else return 0.0;
      }
               break;

      default :
         std::cout << "invalid option. exiting.\n";
         exit(EXIT_FAILURE);

   } /* end of switch */

} /* end of calculate_jij function */

void populate_supercell() {

   /* loop through dimensions */
   for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
         for (int k=0; k<3; k++)

            /* loop through atoms in unitcell */
            for (int atom=0; atom<unitcell.size(); atom++) {

               atom_t temp;
               vec_t uc;
               uc.x = i;
               uc.y = j;
               uc.z = k;

               temp.aid = unitcell[atom].aid;
               temp.element = unitcell[atom].element;
               temp.mat = unitcell[atom].mat;

               /* replicate unitcell atoms */
               temp.pos = unitcell[atom].pos + (uc * ucd);

               /* label unitcell coordinates */
               temp.uc = uc;

               /* dummy values for unneeded struct elements */
               temp.hcat = 0;
               temp.gid = 0;   // this isn't needed yet as the calculations are generated for a large system later

               /* place atom in array */
               supercell.push_back(temp);
            }

   std::cout << "atoms in super cell: "
             << supercell.size() << std::endl;

   array_to_rasmol(supercell, "supercell");

}

