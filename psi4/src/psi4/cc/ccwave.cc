/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "ccwave.h"

#include <vector>
#include <map>

#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

namespace psi {
namespace cc {

CCWavefunction::CCWavefunction(std::shared_ptr<Wavefunction> ref_wfn)
    : Wavefunction(Process::environment.options), params_(Process::environment.options) {
    // Copy the wavefuntion then update
    shallow_copy(ref_wfn);
    set_reference_wavefunction(ref_wfn);
    common_init();
}

CCWavefunction::CCWavefunction(std::shared_ptr<Wavefunction> ref_wfn, Options &options)
    : Wavefunction(options), params_(options) {
    // Copy the wavefuntion then update
    shallow_copy(ref_wfn);
    set_reference_wavefunction(ref_wfn);
    common_init();
}

CCWavefunction::~CCWavefunction() {}

double CCWavefunction::compute_energy() { return 0.0; }

void CCWavefunction::title(std::string &wfn) {
    outfile->Printf("\n");
    outfile->Printf("         ---------------------------------------------------------\n");
    outfile->Printf("                          Coupled Cluster\n");
    outfile->Printf("                           %s wavefunction\n", wfn.c_str());
    outfile->Printf("\n");
    outfile->Printf("                 T. Daniel Crawford\n");
    outfile->Printf("         ---------------------------------------------------------\n");
    outfile->Printf("\n");
}

void CCWavefunction::common_init() {
    title(params_.wfn);
    // get_mo_info();

    // Print out information
    params_.print_parameters(memory_);
}

void CCWavefunction::init_dpd() {
    // cachefiles_.reserve(PSIO_MAXUNIT);
    // std::vector<int *> spaces;
    // std::vector<int *> aospaces;
    // if (params_.ref == 2) {
    //    cachelist_ = cacheprep_uhf(params_.cachelev, cachefiles_.data());
    //    spaces.push_back(moinfo_.aoccpi);
    //    spaces.push_back(moinfo_.aocc_sym);
    //    spaces.push_back(moinfo_.avirtpi);
    //    spaces.push_back(moinfo_.avir_sym);
    //    spaces.push_back(moinfo_.boccpi);
    //    spaces.push_back(moinfo_.bocc_sym);
    //    spaces.push_back(moinfo_.bvirtpi);
    //    spaces.push_back(moinfo_.bvir_sym);
    //    if (params_.aobasis != "NONE") {
    //        aospaces.push_back(moinfo_.aoccpi);
    //        aospaces.push_back(moinfo_.aocc_sym);
    //        aospaces.push_back(moinfo_.sopi);
    //        aospaces.push_back(moinfo_.sosym);
    //        aospaces.push_back(moinfo_.boccpi);
    //        aospaces.push_back(moinfo_.bocc_sym);
    //        aospaces.push_back(moinfo_.sopi);
    //        aospaces.push_back(moinfo_.sosym);
    //    }
    //} else {
    //    cachelist_ = cacheprep_rhf(params_.cachelev, cachefiles_.data());
    //    spaces.push_back(moinfo_.occpi);
    //    spaces.push_back(moinfo_.occ_sym);
    //    spaces.push_back(moinfo_.virtpi);
    //    spaces.push_back(moinfo_.vir_sym);
    //    if (params_.aobasis != "NONE") {
    //        aospaces.push_back(moinfo_.occpi);
    //        aospaces.push_back(moinfo_.occ_sym);
    //        aospaces.push_back(moinfo_.sopi);
    //        aospaces.push_back(moinfo_.sosym);
    //    }
    //}

    // dpd_["mo"].init(0, moinfo_.nirreps, params_.memory, params_.cachetype, cachefiles_.data(), cachelist_,
    //                cache_priority_list_.data(), spaces.size() / 2, spaces);

    // if (aospaces.size()) {
    //    dpd_["ao"].init(1, moinfo_.nirreps, params_.memory, 0, cachefiles_.data(), cachelist_, nullptr,
    //                    aospaces.size() / 2, aospaces);
    //}
}

void CCWavefunction::tear_down() {
    // Free up cache
    //    for (auto &&i : dpd_) {
    //        i.second.file2_cache_close();
    //        i.second.file4_cache_close();
    //    }
    //
    //    if (params_.ref == 2) {
    //        cachedone_uhf(cachelist_);
    //    } else {
    //        cachedone_rhf(cachelist_);
    //    }

    // Clean up I/O
    // exit_io();
}

CCParams::CCParams(Options &options) {
    wfn = options.get_str("WFN");
    if (wfn == "NONE") throw PsiException("Invalid value of input keyword WFN", __FILE__, __LINE__);

    newtrips = options.get_bool("NEW_TRIPLES");

    if (wfn == "BCCD" || wfn == "BCCD_T") {
        brueckner = 1;
    } else {
        brueckner = 0;
    }

    df = (options.get_str("CC_TYPE") == "DF");

    semicanonical = false;
    ref = Reference::RHF;
    auto junk = options.get_str("REFERENCE");
    if (junk == "RHF") {
        ref = Reference::RHF;
    } else if (junk == "ROHF" &&
               (wfn == "MP2" || wfn == "CCSD_T" || wfn == "CCSD_AT" || wfn == "CC3" || wfn == "EOM_CC3" ||
                wfn == "CC2" || wfn == "EOM_CC2" || wfn == "BCCD" || wfn == "BCCD_T")) {
        ref = Reference::UHF;
        semicanonical = true;
    } else if (junk == "ROHF") {
        ref = Reference::ROHF;
    } else if (junk == "UHF") {
        ref = Reference::UHF;
    } else {
        throw PsiException("Invalid value of input keyword REFERENCE", __FILE__, __LINE__);
    }

    // Allow user to force semicanonical
    if (options["SEMICANONICAL"].has_changed()) {
        semicanonical = options.get_bool("SEMICANONICAL");
        ref = Reference::UHF;
    }

    analyze = options.get_bool("ANALYZE");

    dertype = DerivativeType::NONE;
    junk = options.get_str("DERTYPE");
    if (junk == "NONE") {
        dertype = DerivativeType::NONE;
    } else if (junk == "FIRST") {
        dertype = DerivativeType::FIRST;
    } else if (junk == "RESPONSE") {
        dertype = DerivativeType::RESPONSE; /* linear response */
    } else {
        throw PsiException("Invalid value of input keyword DERTYPE", __FILE__, __LINE__);
    }

    print = options.get_int("PRINT");
    maxiter = options.get_int("MAXITER");
    convergence = options.get_double("R_CONVERGENCE");
    e_convergence = options.get_double("E_CONVERGENCE");
    restart = options.get_bool("RESTART");

    aobasis = options.get_str("AO_BASIS");
    cachelevel = options.get_int("CACHELEVEL");

    cachetype = CacheType::LOW;
    junk = options.get_str("CACHETYPE");
    if (junk == "LOW") {
        cachetype = CacheType::LOW;
    } else if (junk == "LRU") {
        cachetype = CacheType::LRU;
    } else {
        throw PsiException("Error in input: invalid CACHETYPE", __FILE__, __LINE__);
    }

    /* No LOW cacheing yet for UHF references */
    if (ref == Reference::UHF) cachetype = CacheType::LRU;

    nthreads = Process::environment.get_n_threads();
    if (options["CC_NUM_THREADS"].has_changed()) {
        nthreads = options.get_int("CC_NUM_THREADS");
    }

    diis = options.get_bool("DIIS");
    t2_coupled = options.get_bool("T2_COUPLED");
    prop = options.get_str("PROPERTY");
    abcd = options.get_str("ABCD");
    local = options.get_bool("LOCAL");
    local_cutoff = options.get_double("LOCAL_CUTOFF");
    local_method = options.get_str("LOCAL_METHOD");
    local_weakp = options.get_str("LOCAL_WEAKP");

    local_cphf_cutoff = options.get_double("LOCAL_CPHF_CUTOFF");
    local_freeze_core = (options.get_str("FREEZE_CORE") != "FALSE");

    local_pairdef = options.get_str("LOCAL_PAIRDEF");
    if (local && dertype == DerivativeType::RESPONSE) {
        local_pairdef = "RESPONSE";
    } else if (local) {
        local_pairdef = "BP";
    }

    num_amps = options.get_int("NUM_AMPS_PRINT");
    bconv = options.get_double("BRUECKNER_ORBS_R_CONVERGENCE");

    // Tying orbital convergence to the desired e_conv,
    //   particularly important for sane numerical frequencies by energy
    if (options["BRUECKNER_ORBS_R_CONVERGENCE"].has_changed()) {
        bconv = options.get_double("BRUECKNER_ORBS_R_CONVERGENCE");
    } else {
        bconv = 100.0 * e_convergence;
    }

    print_mp2_amps = options.get_bool("MP2_AMPS_PRINT");
    print_pair_energies = options.get_bool("PAIR_ENERGIES_PRINT");
    spinadapt_energies = options.get_bool("SPINADAPT_ENERGIES");
    t3_Ws_incore = options.get_bool("T3_WS_INCORE");

    /* get parameters related to SCS-MP2 or SCS-N-MP2 */
    /* see papers by S. Grimme or J. Platz */
    scsn = options.get_bool("SCSN_MP2");
    scs = options.get_bool("SCS_MP2");
    scscc = options.get_bool("SCS_CCSD");
    scsmp2_scale_os = options.get_double("MP2_OS_SCALE");
    scsmp2_scale_ss = options.get_double("MP2_SS_SCALE");
    /* see paper by T. Takatani*/
    scscc_scale_os = options.get_double("CC_OS_SCALE");
    scscc_scale_ss = options.get_double("CC_SS_SCALE");

    if (options["MP2_OS_SCALE"].has_changed() || options["MP2_SS_SCALE"].has_changed()) {
        scs = true;
    }

    if (options["CC_OS_SCALE"].has_changed() || options["CC_SS_SCALE"].has_changed()) {
        scscc = true;
    }
}

void CCParams::print_parameters(size_t memory) const {
    outfile->Printf("\n    Input parameters:\n");
    outfile->Printf("    -----------------\n");
    outfile->Printf("    Wave function   =     %s\n", wfn.c_str());

    if (semicanonical) {
        outfile->Printf("    Reference wfn   =     ROHF changed to UHF for Semicanonical Orbitals\n");
    } else {
        outfile->Printf("    Reference wfn   =     %s\n",
                        (ref == Reference::RHF) ? "RHF" : ((ref == Reference::ROHF) ? "ROHF" : "UHF"));
    }
    outfile->Printf("    Brueckner       =     %s\n", brueckner ? "Yes" : "No");
    if (brueckner) outfile->Printf("    Brueckner conv. =     %3.1e\n", bconv);
    outfile->Printf("    Memory [GiB]    =     %.3f\n", memory * 8 / (1024 * 1024 * 1024.0));
    outfile->Printf("    Maxiter         =    %4d\n", maxiter);
    outfile->Printf("    R_Convergence   =     %3.1e\n", convergence);
    outfile->Printf("    E_Convergence   =     %3.1e\n", e_convergence);
    outfile->Printf("    Restart         =     %s\n", restart ? "Yes" : "No");
    outfile->Printf("    DIIS            =     %s\n", diis ? "Yes" : "No");
    outfile->Printf("    AO Basis        =     %s\n", aobasis.c_str());
    outfile->Printf("    ABCD            =     %s\n", abcd.c_str());
    outfile->Printf("    Cache Level     =     %1d\n", cachelevel);
    outfile->Printf("    Cache Type      =    %4s\n", (cachetype == CacheType::LOW) ? "LOW" : "LRU");
    outfile->Printf("    Print Level     =     %1d\n", print);
    outfile->Printf("    Num. of threads =     %d\n", nthreads);
    outfile->Printf("    # Amps to Print =     %1d\n", num_amps);
    outfile->Printf("    Print MP2 Amps? =     %s\n", print_mp2_amps ? "Yes" : "No");
    outfile->Printf("    Analyze T2 Amps =     %s\n", analyze ? "Yes" : "No");
    outfile->Printf("    Print Pair Ener =     %s\n", print_pair_energies ? "Yes" : "No");

    if (print_pair_energies) outfile->Printf("    Spinadapt Ener. =     %s\n", spinadapt_energies ? "Yes" : "No");
    outfile->Printf("    Local CC        =     %s\n", local ? "Yes" : "No");

    if (wfn == "CC3" || wfn == "EOM_CC3")
        outfile->Printf("    T3 Ws incore    =     %s\n", t3_Ws_incore ? "Yes" : "No");

    if (local) {
        outfile->Printf("    Local Cutoff       =     %3.1e\n", local_cutoff);
        outfile->Printf("    Local Method      =     %s\n", local_method.c_str());
        outfile->Printf("    Weak pairs        =     %s\n", local_weakp.c_str());
        outfile->Printf("    Local pairs       =     %s\n", local_pairdef.c_str());
        outfile->Printf("    Local CPHF cutoff =     %3.1e\n", local_cphf_cutoff);
    }
    outfile->Printf("    SCS-MP2         =     %s\n", scs ? "True" : "False");
    outfile->Printf("    SCSN-MP2        =     %s\n", scsn ? "True" : "False");
    outfile->Printf("    SCS-CCSD        =     %s\n", scscc ? "True" : "False");
    if (scs) {
        outfile->Printf("    SCS_MP2_OS_SCALE =     %.2f\n", scsmp2_scale_os);
        outfile->Printf("    SCS_MP2_SS_SCALE =     %.2f\n", scsmp2_scale_ss);
    }
    if (scsn) {
        outfile->Printf("    SCSN_MP2_OS_SCALE =     %.2f\n", 0.0);
        outfile->Printf("    SCSN_MP2_SS_SCALE =     %.2f\n", 1.76);
    }
    if (scscc) {
        outfile->Printf("    CC_OS_SCALE     =     %.2f\n", scscc_scale_os);
        outfile->Printf("    CC_SS_SCALE     =     %.2f\n", scscc_scale_ss);
    }

    outfile->Printf("\n");
}

}  // namespace cc
}  // namespace psi
