#include <iostream>
#include <fstream>

/* project header files */
#include "./initialise.hpp"
#include "./data.hpp"
#include "./main.hpp"
#include "./io.hpp"
#include "./interactions.hpp"

int calculate_interactions() {

   std::cout << std::endl;
   std::cout << "calculating interactions...\n\n";

   /* first and last index of central unitcell */
   int start;
   int end;

   int interaction_count = 0;
   int tt_interaction_count = 0;
   int rt_interaction_count = 0;

   double total_neighbour_distance = 0;
   double total_tt_neighbour_distance = 0;
   double total_rt_neighbour_distance = 0;

      start = (supercell.size()-unitcell.size())/2;
      end   = (supercell.size()+unitcell.size())/2;

      /* loop through central cell */
      for (int i=start; i<end; ++i) {

         /* loop through super cell (looking for atom j) */
         for (int j=0; j<supercell.size(); ++j) {

            /* create a temporary pair of atoms */
            pair_t pair;
            pair.i = supercell[i];
            pair.j = supercell[j];

            /* if distance less than rcut and not same atom */
            if (pair.are_within_range(sys.rcut_rt, sys.rcut_tt) && !pair.are_same_atom()) {

               /* add neighbour distance to total */
               total_neighbour_distance += pair.rij();

               /* assign interaction id */
               pair.iid = interaction_count;

               /* unitcell displacement */
               pair.ucd = pair.j.uc - pair.i.uc;

               /* calculate exchange energy */
               pair.exchange = calculate_jij(pair);

               /* put interaction into array */
               if (pair.exchange >= 1e-30) {

                  uc_interactions.push_back(pair);
                  interaction_count ++;

                  /* if interaction is fe-fe or fe-co */
                  if (pair.i.is_tm() && pair.j.is_tm()) {
                     tt_interaction_count ++;
                     total_tt_neighbour_distance += pair.rij();
                  }

                  /* if interaction is fe-re */
                  else if ((pair.i.is_re() && pair.j.is_tm()) || (pair.i.is_tm() && pair.j.is_re())) {
                     rt_interaction_count ++;
                     total_rt_neighbour_distance += pair.rij();

                  }
               }
            }
         }
      }

      std::cout
         << "tm-tm interactions: "
         << tt_interaction_count << std::endl;

      std::cout
         << "re-tm interactions: "
         << rt_interaction_count << std::endl;

      std::cout
         << "total interactions: "
         << uc_interactions.size() << std::endl;

      std::cout << std::endl;

      std::cout
         << "mean neighbour distance: "
         << total_neighbour_distance/uc_interactions.size()
         << " A" << std::endl;

      std::cout
         << "mean t-t neighbour distance: "
         << total_tt_neighbour_distance/tt_interaction_count
         << " A" << std::endl;

      std::cout
         << "mean r-t neighbour distance: "
         << total_rt_neighbour_distance/rt_interaction_count
         << " A" << std::endl;

      /* matrix to store interaction strengths */
      std::vector<std::vector<double> > element_interactions;
      element_interactions.resize(elements.size());
      for (int i=0; i<elements.size(); ++i) element_interactions[i].resize(elements.size());

      /* matrix of same dimensions to hold number of interactions between each type */
      std::vector<std::vector<double> > n_interactions;
      n_interactions.resize(elements.size());
      for (int i=0; i<elements.size(); ++i) n_interactions[i].resize(elements.size());

      /* initialise both matrices to zero */
      for (int i=0; i<elements.size(); ++i)
         for (int j=0; j<elements.size(); ++j) {

            element_interactions[i][j] = 0;
            n_interactions[i][j] = 0;

         }

      /* fill matrices */
      for (int interaction=0; interaction<uc_interactions.size(); ++interaction) {

         int i = uc_interactions[interaction].i.mat;
         int j = uc_interactions[interaction].j.mat;

         element_interactions[i][j] += uc_interactions[interaction].exchange;
         n_interactions[i][j] ++;

      }

      /* print matrices */
      std::cout << "\nn_interactions matrix\n";

      for (int i=0; i<elements.size(); ++i) {
         for (int j=0; j<elements.size(); ++j)

            std::cout << n_interactions[i][j] << "\t";

         std::cout << std::endl;
      }

      std::cout << "\nexchange matrix (total)\n\n";

      for (int i=0; i<elements.size(); ++i) {
         for (int j=0; j<elements.size(); ++j) {

            if (n_interactions[i][j] != 0) std::cout << element_interactions[i][j] << "\t";

            else std::cout << "0.000000000" << "\t";

         }

         std::cout << std::endl;
      }

      std::cout << "\nexchange matrix (mean)\n\n";

      for (int i=0; i<elements.size(); ++i) {
         for (int j=0; j<elements.size(); ++j) {
            if (n_interactions[i][j] != 0) {

               element_interactions[i][j] /= n_interactions[i][j];
               std::cout << element_interactions[i][j] << "\t";

            }

            else std::cout << "0.000000000" << "\t";
         }

         std::cout << std::endl;
      }

      std::cout << std::endl;

      outfile
         << "# Interactions n exctype, id i j dx dy dz Jij\n"
         << uc_interactions.size() << "\tisotropic\n";

      for (int i=0; i<uc_interactions.size(); ++i)
         outfile
            << uc_interactions[i].iid       << "\t"
            << uc_interactions[i].i.aid     << "\t"
            << uc_interactions[i].j.aid     << "\t"
            << uc_interactions[i].ucd.x     << "\t"
            << uc_interactions[i].ucd.y     << "\t"
            << uc_interactions[i].ucd.z     << "\t"
            << uc_interactions[i].exchange  << "\n";

      std::cout << "interaction data output to \'output.ucf\'\n";

      return EXIT_SUCCESS;
}

/* calculate exchange energy */
double calculate_jij(pair_t pair) {

   switch (mat.id()) {

      case 1 : {  /* bccfe */
                  /* parameters from Pajda (2001) */
                  /* factor 2 to correct for Hamiltonian */
                  double a = 121.00658;
                  double b = 1.72543313196278;
                  double c = 1e-21;
                  double rij = pair.rij();
                  return c * (a / (rij*rij*rij) - b);
               }
               break;

      case 2 : { /* ndfeb */

                  double a = 121.00658;
                  double b = 1.72543313196278;
                  double c = 1e-21;
                  double rij = pair.rij();

                  double ndfeb_tt_factor = sys.tt_factor;
                  double ndfeb_rt_factor = sys.rt_factor;

                  /* Nd-Nd */
                  if (pair.i.is_re() && pair.j.is_re()) return 0.0;

                  /* Nd-Fe */
                  else if ((pair.i.is_re() && pair.j.is_tm()) || (pair.i.is_tm() && pair.j.is_re()))
                     return sys.rt_constant;

                  /* Fe-Fe */
                  else if (pair.i.is_tm() && pair.j.is_tm())
                     return ndfeb_tt_factor * c*(a/(rij*rij*rij)-b);

                  else return 0.0;

               }
               break;

      case 3 : { // ndfe12

                  double a = 121.00658;
                  double b = 1.72543313196278;
                  double c = 1e-21;
                  double rij = pair.rij();

                  double ndfe12_tt_factor = sys.tt_factor;
                  double ndfe12_rt_factor = sys.rt_factor;

                  /* Nd-Nd */
                  if (pair.i.is_re() && pair.j.is_re()) return 0.0;

                  /* Nd-Fe *** values from matsumoto (2016) *** */
                  else if (((pair.i.is_tm() && pair.j.is_re()) || (pair.i.is_re() && pair.j.is_tm())) && (rij<=sys.rcut_rt))
                     return sys.rt_constant;

                  /* Fe-Fe */
                  else if (pair.i.is_tm() && pair.j.is_tm())
                     return sys.tt_factor * c*(a/(rij*rij*rij)-b);

                  else return 0.0;
               }
               break;

      case 4 : { /* smfe12 */

                  double a = 121.00658;
                  double b = 1.72543313196278;
                  double c = 1e-21;
                  double rij = pair.rij();

                  double smfe12_tt_factor = sys.tt_factor;
                  double smfe12_rt_factor = sys.rt_factor;

                  /* Sm-Sm */
                  if (pair.i.is_re() && pair.j.is_re()) return 0.0;

                  /* Sm-Fe */
                  else if (((pair.i.is_tm() && pair.j.is_re()) || (pair.i.is_re() && pair.j.is_tm())) && (rij<=sys.rcut_rt))
                     return sys.rt_constant;

                  /* Fe-Fe */
                  else if (pair.i.is_tm() && pair.j.is_tm())
                     return sys.tt_factor * c*(a/(rij*rij*rij)-b);

                  else return 0.0;

                  /* these factors come from simulations with re-tm cut-off 4.0 A */
                  // double smfe12_tt_factor = 0.65;
                  // double smfe12_rt_factor = 0.20;

               }
               break;

      case 5 : { /* smco12 */

                  double a = 121.00658;
                  double b = 1.72543313196278;
                  double c = 1e-21;
                  double rij = pair.rij();

                  double smo12_tt_factor = sys.tt_factor;
                  double smo12_rt_factor = sys.rt_factor;

                  /* Sm-Sm */
                  if (pair.i.is_re() && pair.j.is_re()) return 0.0;

                  /* Sm-Co */
                  else if (((pair.i.is_tm() && pair.j.is_re()) || (pair.i.is_re() && pair.j.is_tm())) && (rij<=sys.rcut_rt))
                     return sys.rt_constant;

                  /* Co-Co */
                  else if (pair.i.is_tm() && pair.j.is_tm())
                     return sys.tt_factor * c*(a/(rij*rij*rij)-b);

                  else return 0.0;
               }
               break;

      case 6 : { /* interface */

                  double a = 121.00658;
                  double b = 1.72543313196278;
                  double c = 1e-21;
                  double rij = pair.rij();

                  double ndfeb_tt_factor = sys.tt_factor;
                  double ndfeb_rt_factor = sys.rt_factor;

                  /* Nd-Nd */
                  if (pair.i.is_re() && pair.j.is_re()) return 0.0;

                  /* Nd-Fe */
                  else if ((pair.i.is_re() && pair.j.is_tm()) || (pair.i.is_tm() && pair.j.is_re()))
                     return sys.rt_constant;

                  /* Fe-Fe */
                  else if (pair.i.is_tm() && pair.j.is_tm())
                     return ndfeb_tt_factor * c*(a/(rij*rij*rij)-b);

                  else return 0.0;

               }
               break;

      case 7 : { // ndfeti12

                  double a = 121.00658;
                  double b = 1.72543313196278;
                  double c = 1e-21;
                  double rij = pair.rij();

                  double ndfe12_tt_factor = sys.tt_factor;
                  double ndfe12_rt_factor = sys.rt_factor;

                  /* Nd-Nd */
                  if (pair.i.is_re() && pair.j.is_re()) return 0.0;

                  /* Nd-Fe *** values from matsumoto (2016) *** */
                  else if (((pair.i.is_tm() && pair.j.is_re()) || (pair.i.is_re() && pair.j.is_tm())) && (rij<=sys.rcut_rt))
                     return sys.rt_constant;

                  /* Fe-Fe */
                  else if (pair.i.is_tm() && pair.j.is_tm())
                     return sys.tt_factor * c*(a/(rij*rij*rij)-b);

                  else return 0.0;
               }
               break;

      default :
               std::cout << "no exchange function found for this material. exiting.\n";
               exit(EXIT_FAILURE);

   } /* end of switch */

} /* end of calculate_jij function */

int populate_supercell() {

   std::cout << "\ncreating and populating super cell...\n";

   int global_id_counter = 0;

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
               temp.pos = unitcell[atom].pos + (uc * mat.ucd);

               /* label unitcell coordinates */
               temp.uc = uc;

               /* dummy values for unneeded struct elements */
               temp.hcat = 0;
               temp.gid = global_id_counter;
               global_id_counter ++;

               /* place atom in array */
               supercell.push_back(temp);
            }

   std::cout
      << "atoms in super cell: " << supercell.size() << "\n\n";

//   array_to_rasmol(supercell, "supercell");

   return EXIT_SUCCESS;
}

