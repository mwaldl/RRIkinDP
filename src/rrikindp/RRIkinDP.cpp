/**
 *
 * RRIkinDP: Model RNA RNA interaction kinetics on direct paths.
 *
 * Copyright (C) Maria Waldl <maria@waldl.org> <maria@tbi.univie.ac.at>
 *
 */

#include <iostream>
#include <string>
// #include <stdlib>
// #include <cstring>
#include <fstream>
#include <limits>
#include <vector>
#include <cstdlib> // import free()

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/program_options.hpp>
#include <easylogging++.h>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/spirit/home/x3.hpp>

#include <IntaRNA/Accessibility.h>
#include <IntaRNA/AccessibilityVrna.h>
#include <IntaRNA/Interaction.h>
#include <IntaRNA/RnaSequence.h>

// #include "IntaRNA/ReverseAccessibility.h"
#include <IntaRNA/AccessibilityConstraint.h>
#include <IntaRNA/InteractionRange.h>
#include <IntaRNA/VrnaHandler.h>

#include <IntaRNA/InteractionEnergy.h>
#include <IntaRNA/InteractionEnergyVrna.h>

extern "C" {
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/MEA.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/constraints/hard.h>
}

INITIALIZE_EASYLOGGINGPP

const double d_infinity = std::numeric_limits<int>::max();

int
map_2D_linear(int k, int l, int matrix_dimension) {
    return k * matrix_dimension + l;
}

double
smallest(double x, double y, double z) {
    return std::min(std::min(x, y), z);
}

class EM {
public:
    EM(const IntaRNA::Interaction,
       const std::string,
       const std::string,
       std::string,
       std::string,
       bool,
       bool,
       double);
    double
    get_e(int, int) const;
    double
    get_hybride_e(int, int) const;
    double
    get_accessibility(int, int) const;
    double
    get_ED1(int, int) const;
    double
    get_ED2(int, int) const;
    std::string
    get_structure1() const;
    std::string
    get_structure2() const;
    std::string
    get_structures() const;

private:
    IntaRNA::Interaction interaction_;
    int length_I_;
    //   void set_e (int, int, double);
    std::vector<double> matrix_;
    std::vector<double> loop_e_;
    std::vector<double> ext_left_e_;
    std::vector<double> ext_right_e_;
    std::vector<double> dangle_left_e_;
    std::vector<double> dangle_right_e_;
    std::vector<double> end_e_;
    std::vector<double> hybrid_e_;
    std::vector<double> accessibility1_e_;
    std::vector<double> accessibility2_e_;
    std::vector<double> accessibility_e_;
    void
    set_loop_e();
    void
    set_ext_loops_e();
    void
    set_hybrid_e();
    void
    set_accessibility_pf();
    void
    set_accessibility_fixed_struct();
    IntaRNA::RnaSequence s1_;
    std::string str1_;
    std::string constraint_str1_;
    IntaRNA::RnaSequence s2_;
    std::string str2_;
    std::string constraint_str2_;
    double temperature_;
    IntaRNA::VrnaHandler handler_;
    IntaRNA::AccessibilityVrna access1_;
    IntaRNA::AccessibilityVrna access2_;
    IntaRNA::ReverseAccessibility rev_access2_;
    IntaRNA::InteractionEnergyVrna energy_;
    bool dangle_;
    bool use_pf_;
};

EM::EM(IntaRNA::Interaction input_interaction,
       std::string s_a,
       std::string s_b,
       std::string str_a = "",
       std::string str_b = "",
       bool accessibilty_from_pf = true,
       bool dangle = true,
       double temperature = 37.0)
    : interaction_(input_interaction),
      length_I_(interaction_.basePairs.size()),
      matrix_(length_I_ * length_I_, d_infinity),
      loop_e_(length_I_ - 1, d_infinity),
      ext_left_e_(length_I_, d_infinity),
      ext_right_e_(length_I_, d_infinity),
      dangle_left_e_(length_I_, d_infinity),
      dangle_right_e_(length_I_, d_infinity),
      end_e_(length_I_, d_infinity),
      hybrid_e_(length_I_ * length_I_, d_infinity),
      accessibility1_e_(length_I_ * length_I_, d_infinity),
      accessibility2_e_(length_I_ * length_I_, d_infinity),
      accessibility_e_(length_I_ * length_I_, d_infinity),
      s1_("A", s_a),
      str1_(str_a),
      constraint_str1_(""),
      s2_("B", s_b),
      str2_(str_b),
      constraint_str2_(""),
      temperature_(temperature),
      handler_(temperature_, "Turner04", false, false), // TODO
      access1_(s1_, s1_.size(), NULL, handler_, s1_.size()),
      access2_(s2_, s2_.size(), NULL, handler_, s2_.size()),
      rev_access2_(access2_),
      energy_(access1_, rev_access2_, handler_, 30, 30, false),
      use_pf_(accessibilty_from_pf),
      dangle_(dangle) {
    set_loop_e();
    set_ext_loops_e();
    set_hybrid_e();
    if (use_pf_) {
        set_accessibility_pf();
    } else {
        set_accessibility_fixed_struct();
    }
}

void
EM::set_loop_e() {
    for (int i = 0; i < length_I_ - 1; i++) {
        int i1 = interaction_.basePairs[i].first;
        int j1 = interaction_.basePairs[i + 1].first;
        int i2 = interaction_.basePairs[i].second;
        int j2 = interaction_.basePairs[i + 1].second;
        int j2_rev = rev_access2_.getReversedIndex(j2);
        int i2_rev = rev_access2_.getReversedIndex(i2);
        auto x = energy_.getE_interLeft(i1, j1, i2_rev, j2_rev);
        //  i1 interacts with i2, j1 interacts with j2; 1 == first sequence, 2
        //  == second sequence
        loop_e_[i] = x;
    }
}

void
EM::set_ext_loops_e() {
    for (int i = 0; i < length_I_; i++) {
        int i1 = interaction_.basePairs[i].first;
        int i2 = interaction_.basePairs[i].second;
        int i2_rev = rev_access2_.getReversedIndex(i2);
        int dangleL_e = 0;
        int dangleR_e = 0;
        if (dangle_) {
            dangleL_e = energy_.getE_danglingLeft(i1, i2_rev);
            dangleR_e = energy_.getE_danglingRight(i1, i2_rev);
        }
        int end_e = energy_.getE_endLeft(i1, i2_rev);
        ext_left_e_[i] = dangleL_e + end_e;
        ext_right_e_[i] = dangleR_e + end_e;
        dangle_left_e_[i] = dangleL_e;
        dangle_right_e_[i] = dangleR_e;
        end_e_[i] = end_e;
    }
}

void
EM::set_hybrid_e() {
    double init_e = energy_.getE_init();

    for (int k = 0; k < length_I_; k++) {
        for (int l = k; l < length_I_; l++) {
            double h_energy = 0.0;
            for (int i = k; i < l; i++) {
                h_energy += loop_e_[i];
            }
            hybrid_e_[map_2D_linear(k, l, length_I_)] =
                h_energy + ext_left_e_[k] + ext_right_e_[l] + init_e;
            /*
            // IntaRNA like energy
            // ===================

            hybrid_e_[map_2D_linear(k, l, length_I_)] = h_energy +
                dangle_left_e_[k] *
                    energy_.getPr_danglingLeft(i1, j1, i2_rev, j2_rev) +
                end_e_[k] +
                dangle_right_e_[l] *
                    energy_.getPr_danglingRight(i1, j1, i2_rev, j2_rev) +
                end_e_[l] + init_e;
            */
        }
    }
}

void
EM::set_accessibility_pf() {
    // std::cout << "i1, j1, a1, i2, j2, a2"  << std::endl;
    for (int k = 0; k < length_I_; k++) {
        int i1 = interaction_.basePairs[k].first;
        int i2 = interaction_.basePairs[k].second;

        for (int l = k; l < length_I_; l++) {
            int j1 = interaction_.basePairs[l].first;
            int j2 = interaction_.basePairs[l].second;
            auto a1 = access1_.getED(i1, j1);
            auto a2 = access2_.getED(j2, i2);
            double energy = a1 + a2;
            accessibility_e_[map_2D_linear(k, l, length_I_)] = energy;
            accessibility1_e_[map_2D_linear(k, l, length_I_)] = a1;
            accessibility2_e_[map_2D_linear(k, l, length_I_)] = a2;
        }
    }
}

void
EM::set_accessibility_fixed_struct() {
    // set up fold compounds
    // =====================
    vrna_md_t md; // model details
    vrna_md_set_default(&md);
    vrna_fold_compound_t *vc1;
    vc1 = vrna_fold_compound(s1_.asString().c_str(), &md, VRNA_OPTION_MFE);
    vrna_fold_compound_t *vc2;
    vc2 = vrna_fold_compound(s2_.asString().c_str(), &md, VRNA_OPTION_MFE);

    // get initial intramolecular structures as pair tables
    // ====================================================
    short *pt1; // pair table of first RNA
    int e1;     // energy of start structure of first RNA
    short *pt2; // pair table of second RNA
    int e2;     // energy of start structure of second RNA

    if (str1_.empty()) { // use mfe struture
        char *structure = (char *)vrna_alloc(sizeof(char) * (s1_.size() + 1));
        e1 = vrna_mfe(vc1, structure) * 100; // x100 to convert to deca-calories
        pt1 = vrna_ptable(structure);
        str1_ = structure;
    } else { // use input structure
        const char *structure = str1_.data();
        pt1 = vrna_ptable(structure);
        // e1 = vrna_eval_structure_pt(vc1, pt1);
    }

    if (str2_.empty()) { // use mfe struture
        char *structure = (char *)vrna_alloc(sizeof(char) * (s2_.size() + 1));
        e2 = vrna_mfe(vc2, structure) * 100; // x100 to convert to deca-calories
        pt2 = vrna_ptable(structure);
        str2_ = structure;

    } else { // use input structure
        const char *structure = str2_.data();
        pt2 = vrna_ptable(structure);
        // e2 = vrna_eval_structure_pt(vc2, pt2);
    }

    // get unpairing costs for all interaction ranges
    // ==============================================
    for (int k = 0; k < length_I_; k++) {
        short *pt1_refolded = (short *)vrna_alloc(sizeof(short) * (pt1[0] + 1));
        std::copy(pt1, pt1 + pt1[0] + 1, pt1_refolded);
        int cost1 = 0;
        short *pt2_refolded = (short *)vrna_alloc(sizeof(short) * (pt2[0] + 1));
        std::copy(pt2, pt2 + pt2[0] + 1, pt2_refolded);
        int cost2 = 0;
        // continue moving to zero based bp list here
        int i1 = interaction_.basePairs[k].first + 1;
        int i2 = interaction_.basePairs[k].second + 1;
        int old_j1 = i1 - 1;
        int old_j2 = i2 + 1;

        for (int l = k; l < length_I_; l++) {
            int j1 = interaction_.basePairs[l].first + 1;
            int j2 = interaction_.basePairs[l].second + 1;
            for (int j = old_j1 + 1; j <= j1; j++) {
                if (pt1_refolded[j] != 0) {
                    int cost =
                        vrna_eval_move_pt(vc1, pt1_refolded,
                                          std::max(-j, -pt1_refolded[j]),
                                          std::min(-j, -pt1_refolded[j]));
                    cost1 += cost;
                    pt1_refolded[pt1_refolded[j]] = 0;
                    pt1_refolded[j] = 0;
                }
            }
            for (int j = old_j2 - 1; j >= j2; j--) {
                if (pt2_refolded[j] != 0) {
                    int cost =
                        vrna_eval_move_pt(vc2, pt2_refolded,
                                          std::max(-j, -pt2_refolded[j]),
                                          std::min(-j, -pt2_refolded[j]));
                    cost2 += cost;
                    pt2_refolded[pt2_refolded[j]] = 0;
                    pt2_refolded[j] = 0;
                }
            }
            old_j1 = j1;
            old_j2 = j2;
            double energy = cost1 + cost2;
            accessibility_e_[map_2D_linear(k, l, length_I_)] = energy;
            accessibility1_e_[map_2D_linear(k, l, length_I_)] = cost1;
            accessibility2_e_[map_2D_linear(k, l, length_I_)] = cost2;
        }
    }
}

double
EM::get_e(int k, int l) const {
    return hybrid_e_[map_2D_linear(k, l, length_I_)] +
        accessibility_e_[map_2D_linear(k, l, length_I_)];
}

double
EM::get_hybride_e(int k, int l) const {
    int index_1D = map_2D_linear(k, l, length_I_);
    return hybrid_e_[index_1D];
}

double
EM::get_accessibility(int k, int l) const {
    int index_1D = map_2D_linear(k, l, length_I_);
    return accessibility_e_[index_1D];
}

double
EM::get_ED1(int k, int l) const {
    int index_1D = map_2D_linear(k, l, length_I_);
    return accessibility1_e_[index_1D];
}

double
EM::get_ED2(int k, int l) const {
    int index_1D = map_2D_linear(k, l, length_I_);
    return accessibility2_e_[index_1D];
}

std::string
EM::get_structure1() const {
    return str1_;
}
std::string
EM::get_structure2() const {
    return str2_;
}

std::string
EM::get_structures() const {
    std::string faster_header =
        ">" + interaction_.s1->getId() + "&" + interaction_.s2->getId();

    std::string sequences =
        interaction_.s1->asString() + "&" + interaction_.s2->asString();

    // get interaction ranges
    // ======================
    int interaction_start_1 = interaction_.basePairs[0].first;
    int interaction_end_1 = interaction_.basePairs[length_I_ - 1].first;
    int interaction_start_2 = interaction_.basePairs[length_I_ - 1].second;
    int interaction_end_2 = interaction_.basePairs[0].second;

    // inter interaction structure
    // ===========================
    std::string interaction1(interaction_.s1->size(), '.');
    std::string interaction2(interaction_.s2->size(), '.');
    for (int i = 0; i < length_I_; i++) {
        interaction1[interaction_.basePairs[i].first] = '(';
        interaction2[interaction_.basePairs[i].second] = ')';
        // std::cout << i << " " << interaction_.basePairs[i].first << " "
        //          << interaction_.basePairs[i].second << std::endl;
    }
    std::string inter_str = interaction1 + "&" + interaction2;

    // set up fold compounds
    // =====================
    vrna_md_t md; // model details
    vrna_md_set_default(&md);
    vrna_fold_compound_t *vc1;
    vc1 = vrna_fold_compound(interaction_.s1->asString().c_str(), &md,
                             VRNA_OPTION_DEFAULT);
    vrna_fold_compound_t *vc2;
    vc2 = vrna_fold_compound(interaction_.s2->asString().c_str(), &md,
                             VRNA_OPTION_DEFAULT);

    // get intramolecular structures
    // =============================
    if (use_pf_) {
        // get mfe
        // =======

        double mfe1;
        double mfe2;

        char *structure1 = (char *)vrna_alloc(sizeof(char) * (s1_.size() + 1));
        mfe1 =
            vrna_mfe(vc1,
                     structure1); // *100; // x100 to convert to deca-calories
        std::string intra_str1_mfe = structure1;
        free(structure1);
        char *structure2 = (char *)vrna_alloc(sizeof(char) * (s2_.size() + 1));
        mfe2 =
            vrna_mfe(vc2,
                     structure2); // *100; // x100 to convert to deca-calories
        std::string intra_str2_mfe = structure2;
        free(structure2);
        std::string start_intra_str = intra_str1_mfe + "&" +   intra_str2_mfe;


        // get initial intramolecular mea structures
        // =========================================
        vrna_exp_params_rescale(vc1, &mfe1);
        vrna_exp_params_rescale(vc2, &mfe2);

        char *prob_string1 =
            (char *)vrna_alloc(sizeof(char) * (s1_.size() + 1));
        float en1 = vrna_pf(vc1, prob_string1);
        char *prob_string2 =
            (char *)vrna_alloc(sizeof(char) * (s2_.size() + 1));
        float en2 = vrna_pf(vc2, prob_string2);
        std::string intra_str1_prob = prob_string1;
        std::string intra_str2_prob = prob_string2;
        free(prob_string1);
        free(prob_string2);
        std::string intra_str_prob_start = intra_str1_prob + '&' + intra_str2_prob;

        float mea1;
        char *mea_structure1 = vrna_MEA(vc1, 0.5, &mea1);
        float mea2;
        char *mea_structure2 = vrna_MEA(vc2, 0.5, &mea2);
        //std::cout << mea_structure1 << "&" << mea_structure2 << std::endl;
        std::string intra_str1_mea = mea_structure1;
        std::string intra_str2_mea = mea_structure2;
        free(mea_structure1);
        free(mea_structure2);
        std::string intra_str_mea_start = intra_str1_mea + '&' + intra_str2_mea;

        // constraint strings
        // ======================
        std::string constraint1(interaction_.s1->size(), '.');
        std::string constraint2(interaction_.s2->size(), '.');
        for (int i = interaction_start_1; i <= interaction_end_1; i++) {
            constraint1[i] = 'x';
        }
        for (int i = interaction_start_2; i <= interaction_end_2; i++) {
            constraint2[i] = 'x';
        }
        //std::cout << constraint1 << "&" << constraint2 << std::endl;

        // get constraint intramolecular mea structure
        // ===========================================
        // const char *  	constraint,
        vrna_hc_add_from_db(vc1, constraint1.c_str(),
                            VRNA_CONSTRAINT_DB_DEFAULT);

        char *xprob_string1 =
            (char *)vrna_alloc(sizeof(char) * (s1_.size() + 1));
        en1 = vrna_pf(vc1, xprob_string1);

        vrna_hc_add_from_db(vc2, constraint2.c_str(),
                            VRNA_CONSTRAINT_DB_DEFAULT);

        char *xprob_string2 =
            (char *)vrna_alloc(sizeof(char) * (s2_.size() + 1));
        en2 = vrna_pf(vc2, xprob_string2);

        std::string intra_str1_prob_full = xprob_string1;
        std::string intra_str2_prob_full = xprob_string2;
        free(xprob_string1);
        free(xprob_string2);
        std::string intra_str_prob_full = intra_str1_prob_full + '&' + intra_str2_prob_full;


        float xmea1;
        char *xmea_structure1 = vrna_MEA(vc1, 0.5, &xmea1);
        float xmea2;
        char *xmea_structure2 = vrna_MEA(vc2, 0.5, &xmea2);

        std::string intra_str1_mea_full = xmea_structure1;
        std::string intra_str2_mea_full = xmea_structure2;

        free(xmea_structure1);
        free(xmea_structure2);

        std::string intra_str_mea_full = intra_str1_mea_full + '&' + intra_str2_mea_full;

        // get constraint intamolecular mfe structure
        // ==========================================
        double mfe_const1;
        char *structure1_mfe = (char *)vrna_alloc(sizeof(char) * (s1_.size() + 1));
        mfe_const1 =
            vrna_mfe(vc1,
                     structure1_mfe); // *100; // x100 to convert to deca-calories
        std::string intra_str1_mfe_full = structure1_mfe;
        free(structure1_mfe);

        double mfe_const2;
        char *structure2_mfe = (char *)vrna_alloc(sizeof(char) * (s2_.size() + 1));
        mfe_const2 =
            vrna_mfe(vc2,
                     structure2_mfe); // *100; // x100 to convert to deca-calories
        std::string intra_str2_mfe_full = structure2_mfe;
        free(structure2_mfe);

        std::string intra_str_mfe_full = intra_str1_mfe_full + '&' + intra_str2_mfe_full;

        // return fasta string
        // ===================
        return faster_header + " start mfe\n" +
            sequences + "\n" +
            start_intra_str +"\n"  +
            faster_header + " start mea\n" +
            sequences + "\n" +
            intra_str_mea_start +"\n"  +
            faster_header + " start prob\n" +
            sequences + "\n" +
            intra_str_prob_start + "\n" +
            faster_header + " full mfe\n" +
            sequences + "\n" +
            intra_str_mfe_full + "\n" +
            inter_str + "\n" +
            faster_header + " full mea\n" +
            sequences + "\n" +
            intra_str_mea_full + "\n" +
            inter_str + "\n" +
            faster_header + " full prob\n" +
            sequences + "\n" +
            intra_str_prob_full + "\n" +
            inter_str + "\n";

    } else {
        const char *structure1 = str1_.data();
        short *pt1 = vrna_ptable(structure1);
        for (int i = interaction_start_1 + 1; i <= interaction_end_1 + 1; i++) {
            if (pt1[i] != 0) {
                pt1[pt1[i]] = 0;
                pt1[i] = 0;
            }
        }
        const char *structure2 = str2_.data();
        short *pt2 = vrna_ptable(structure2);
        for (int i = interaction_start_2 + 1; i <= interaction_end_2 + 1; i++) {
            if (pt2[i] != 0) {
                pt2[pt2[i]] = 0;
                pt2[i] = 0;
                /*
                #include <ViennaRNA/landscape/move.h>
                            void vrna_move_apply 	( 	short *
                pt, const vrna_move_t *  	m
                        )

                vrna_move_t *  (pos5, pos3, next(pointer to next move))

                  */
            }
        }
        std::string intra_str1 = vrna_db_from_ptable(pt1);
        std::string intra_str2 = vrna_db_from_ptable(pt2);
        std::string intra_str = intra_str1 + "&" + intra_str2;

        std::string start_intra_str = str1_ + "&" + str2_;

        // return fasta string
        // ===================
        return faster_header + " start\n" + sequences + "\n" + start_intra_str +
            "\n" + faster_header + " full\n" + sequences + "\n" + intra_str +
            "\n" + inter_str + "\n";
    }
}

double
minBarrier(int start, int len, int seed_len, const EM &em) {
    double barrier[len + 1][len + 1];
    for (int i = start; i >= 0; i--) {
        for (int j = start + seed_len - 1; j < len; j++) {
            //       std::cout << "i: " << i  << ", j: " << j << std::endl;
            double x = d_infinity;
            double y = d_infinity;
            double z = d_infinity;
            if ((start + seed_len - 1 == j) && (start == i)) {
                x = em.get_e(i, j);
            }
            if (i < start) {
                y = std::max(em.get_e(i, j), barrier[i + 1][j]);
            }
            if (j > start + seed_len - 1) {
                z = std::max(em.get_e(i, j), barrier[i][j - 1]);
            }
            barrier[i][j] = smallest(x, y, z);
            //        std::cout  << "barrier: "<< barrier[i][j] << ", energy:
            //        "<< test[i*len+j] << ", b(i+1,j): "<<  barrier[i+1][j] <<
            //        ", b(i,j-1): "<< barrier[i][j-1] << ", x: "<< x << ", y:
            //        "<<  y <<
            //        ", z: "<< z << std::endl;
        }
    }
    return barrier[0][len - 1];
}

//=====================================================================================

int
main(int argc, char **argv) {
    // set overall logging style
    el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Format,
                                       std::string("# %level : %msg"));
    // no log file output
    el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToFile,
                                       std::string("false"));
    // set additional logging flags
    el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
    el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
    el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
    el::Loggers::addFlag(el::LoggingFlag::AllowVerboseIfModuleNotSpecified);

    // setup logging with given parameters
    START_EASYLOGGINGPP(argc, argv);

    namespace po = boost::program_options;
    namespace x3 = boost::spirit::x3;

    // read input
    // =============================

    po::options_description description("DirectPaths Usage");

    description.add_options()("help", "Display this help message")(
        "version", "Display the version number")(
        "id_a", po::value<std::string>()->required(),
        "id of first sequence")("seq_a", po::value<std::string>()->required(),
                                "first sequence")(
        "str_a", po::value<std::string>()->default_value(""),
        "intramolecular structure of first sequence in dotbraket notation")(
        "id_b", po::value<std::string>()->required(),
        "id of second sequence")("seq_b", po::value<std::string>()->required(),
                                 "second sequence")(
        "str_b", po::value<std::string>()->default_value(""),
        "intramolecular structure of second sequence in dotbraket notation")(
        "interaction_bps", po::value<std::string>()->required(),
        "interaction base pair list as string (one based)")(
        "seed", po::value<int>()->required(),
        "seed length")("write_all_barriers", po::value<std::string>(),
                       "file paths to write minmal barriers for all seeds to")(
        "write_states", po::value<std::string>(),
        "file paths to write states and their energies to")(
        "compute_states_only",
        "compute states and output their energies to path "
        "specified in 'write_states'")(
        "no_dangle", "turn off dangle contributions at interaction ends")(
        "fixed_intramolecular_structures",
        "compute accessibilities based on fixed intramolecular structures "
        "instead of based on partition function")(
        "write_structures", po::value<std::string>(),
        "file paths to write intramolecular and fully extended intermolecular "
        "structures to")("temperature",
                         po::value<double>()->default_value(37.0),
                         "temperature in Celsius");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(description).run(),
              vm);

    if (vm.count("help")) {
        std::cout << description << std::endl;
        return 0;
    }


    if (vm.count("version")) {
        std::cout << "v0.0.1" << std::endl;
        return 0;
    }

    po::notify(vm);

    std::string seq_a = vm["seq_a"].as<std::string>();
    std::string seq_a_id = vm["id_a"].as<std::string>();
    std::string seq_b = vm["seq_b"].as<std::string>();
    std::string seq_b_id = vm["id_b"].as<std::string>();
    // std::cout << "before first dtructure" << std::endl;
    std::string bp_string = vm["interaction_bps"].as<std::string>();
    int seed_len = vm["seed"].as<int>();
    double temp = vm["temperature"].as<double>();
    std::string barriers_file;
    bool write_barriers;
    if (vm.count("write_all_barriers")) {
        barriers_file = vm["write_all_barriers"].as<std::string>();
        write_barriers = true;
    } else {
        write_barriers = false;
    }

    bool use_pf = true;
    std::string str_a = vm["str_a"].as<std::string>();
    std::string str_b = vm["str_b"].as<std::string>();
    if (vm.count("fixed_intramolecular_structures")) {
        use_pf = false;
    }
    bool use_dangle = true;
    if (vm.count("no_dangle")) {
        use_dangle = false;
    }

    bool write_states;
    std::string states_file;
    if (vm.count("write_states")) {
        states_file = vm["write_states"].as<std::string>();
        write_states = true;
    } else {
        write_states = false;
    }
    bool write_states_only = false;
    if (vm.count("compute_states_only") && !write_states) {
        std::cout << "no output file path for states provided" << std::endl;
        return 1;
    } else if (vm.count("compute_states_only")) {
        write_states_only = true;
    }

    bool write_structures;
    std::string structures_file;
    if (vm.count("write_structures")) {
        structures_file = vm["write_structures"].as<std::string>();
        write_structures = true;
    } else {
        write_structures = false;
    }

    // process bp_list to get Intarna pairing vector
    //===============================================

    const auto pair = x3::lit('(') > x3::ulong_ > ',' > x3::ulong_ > ')';
    // const auto list = x3::lit('[') > (pair % ',') > ']';
    const auto list = (pair % ':');
    auto iter = bp_string.begin();
    auto end_iter = bp_string.end();
    std::vector<std::pair<long unsigned int, long unsigned int>> bps_list_b1;
    x3::parse(iter, end_iter, list, bps_list_b1);
    std::vector<std::pair<long unsigned int, long unsigned int>> bps_list_b0;

    std::copy(bps_list_b1.begin(), bps_list_b1.end(),
              back_inserter(bps_list_b0));
    for (int i = 0; i < bps_list_b0.size(); i++) {
        bps_list_b0[i].first -= 1;
        bps_list_b0[i].second -= 1;
    }

    // set up interaction
    //====================

    IntaRNA::RnaSequence s1(seq_a_id, seq_a);
    IntaRNA::RnaSequence s2(seq_b_id, seq_b);
    IntaRNA::Interaction interaction(s1, s2);
    interaction.basePairs = bps_list_b0;

    // check input
    //=============
    if (interaction.basePairs.size() < seed_len) {
        std::cout << "seed longer than interaction" << std::endl;
        //       throw std::exception();
        return 1;
    }
    if (!interaction.isValid()) {
        std::cout << "interaction not valid" << std::endl;
        return 1;
    }

    // fill energy matrix - initiate interaction class
    // ==============================================

    EM em(interaction, seq_a, seq_b, str_a, str_b, use_pf, use_dangle, temp);

    int q = interaction.basePairs.size(); // number of base pairs in structure
    double full_energy = em.get_e(0, q - 1);
    double full_hybrid_energy = em.get_hybride_e(0, q - 1);
    int v = q * q + q; // size of energy matrix/array


    // output structures
    // =================

    if (write_structures) {
        std::ofstream structures_file_handler(structures_file);
        structures_file_handler << em.get_structures();
    }
    // TODO: use fasta class instead



    // output states and their energies
    // ================================

    if (write_states) {
        std::ofstream states_file_handler(states_file);
        states_file_handler << "k\tl\tE\tEhybrid\tED1\tED2\tED" << std::endl;
        for (int k = 0; k < q; k++) {
            for (int l = k; l < q; l++) {
                states_file_handler << k << "\t" << l << "\t" << em.get_e(k, l)
                                    << "\t" << em.get_hybride_e(k, l) << "\t"
                                    << em.get_ED1(k, l) << "\t"
                                    << em.get_ED2(k, l) << "\t"
                                    << em.get_accessibility(k, l) << std::endl;
            }
        }
    }
    if (write_states_only) {
        return 0;
    }



    // set up barrier tracking and output
    // ==================================

    int best_start_bp = 0;
    double barrier = d_infinity;
    if (write_barriers) {
        std::ofstream barriers_file_handler(barriers_file);
        barriers_file_handler
            << "start\tmin_barrier\tseed_E\tseed_Ehybrid\tseed_"
               "ED\tseed_ED1\tseed_ED2\tfull_E"
            << std::endl;
    }

    // iterate over all start points (seeds)
    // =====================================

    for (int start = 0; start <= interaction.basePairs.size() - seed_len;
         start++) {
        //   compute min barrier for current seed
        //   ====================================

        double min_barrier = minBarrier(start, q, seed_len, em);

        //   track min barrier
        //   =================

        if (min_barrier < barrier) {
            best_start_bp = start;
            barrier = min_barrier;
        }

        //   write  barrier energy for current seed to barrier file
        //   ======================================================
        if (write_barriers) {
            std::ofstream barriers_file_handler(barriers_file, std::ios::app);
            barriers_file_handler << start << "\t" << min_barrier << "\t";
            barriers_file_handler << em.get_e(start, start + seed_len - 1)
                                  << "\t";
            barriers_file_handler
                << em.get_hybride_e(start, start + seed_len - 1) << "\t";
            barriers_file_handler
                << em.get_accessibility(start, start + seed_len - 1);
            barriers_file_handler
                << "\t" << em.get_ED1(start, start + seed_len - 1) << "\t";
            barriers_file_handler << em.get_ED2(start, start + seed_len - 1)
                                  << "\t";
            barriers_file_handler << full_energy << std::endl;
        }
    }

    // output summary to stdout
    // ========================

    if (use_pf != true) {
        std::cout << "|: id a: " << seq_a_id << std::endl;
        std::cout << "|: seq a: " << seq_a << std::endl;
        std::cout << "|: str a: " << em.get_structure1() << std::endl;
        std::cout << "|: id b: " << seq_b_id << std::endl;
        std::cout << "|: seq b: " << seq_b << std::endl;
        std::cout << "|: str b: " << em.get_structure2() << std::endl;
    }
    std::cout << "|: interaction bps: " << bp_string << std::endl;
    std::cout << "|: interaction size: " << q << std::endl;
    std::cout << "|: full interaction energy: " << full_energy << std::endl;
    std::cout << "|: hybridisation energy: " << full_hybrid_energy << std::endl;
    std::cout << "|: ED1+ED2: " << em.get_accessibility(0, q - 1) << std::endl;
    std::cout << "|: ED1: " << em.get_ED1(0, q - 1) << std::endl;
    std::cout << "|: ED2: " << em.get_ED2(0, q - 1) << std::endl;
    std::cout << "|: best start: " << best_start_bp << std::endl;
    std::cout << "|: lowest barrier: " << barrier << std::endl;
    std::cout << "|: seed length: " << seed_len << std::endl;
    std::cout << "|: seed ED of lowest barrier seed: "
              << em.get_accessibility(best_start_bp,
                                      best_start_bp + seed_len - 1)
              << std::endl;
    std::cout << "|: seed hybridisation energy of lowest barrier seed: "
              << em.get_hybride_e(best_start_bp, best_start_bp + seed_len - 1)
              << std::endl;
    std::cout << "|: seed energy of lowest barrier seed: "
              << em.get_e(best_start_bp, best_start_bp + seed_len - 1)
              << std::endl;

    return 0;
}
