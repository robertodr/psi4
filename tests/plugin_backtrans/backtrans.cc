#include "psi4-dec.h"
#include <libmints/mints.h>
#include <liboptions/liboptions.h>
#include "libchkpt/chkpt.h"
#include "libiwl/iwl.h"
#include "libmints/mints.h"
#include "libtrans/integraltransform.h"
#include "libiwl/iwl.hpp"
#include "libplugin/plugin.h"
#include "ccfiles.h"

#define ID(x) ints.DPD_ID(x)
INIT_PLUGIN

using namespace boost;

namespace psi{ namespace tpdmtest{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "PLUGIN_BACKTRANS"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" PsiReturnType
plugin_backtrans(Options &options)
{
    int print = options.get_int("PRINT");

    shared_ptr<PSIO> psio(_default_psio_lib_);


    std::vector<shared_ptr<MOSpace> > spaces;
    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));
  
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(chkpt, spaces, IntegralTransform::Restricted,
             IntegralTransform::DPDOnly, IntegralTransform::QTOrder, IntegralTransform::None, false);
    dpd_set_default(ints.get_dpd_id());
    ints.set_print(print);
    ints.set_keep_dpd_so_ints(1);
    ints.set_keep_iwl_so_ints(1);
    ints.initialize();
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    ints.backtransform_density();
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    psio->open(PSIF_TPDM_PRESORT, PSIO_OPEN_OLD);
    psio->open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);
    psio->open(PSIF_TPDM_HALFTRANS, PSIO_OPEN_OLD);


    char **labels = Process::environment.molecule()->irrep_labels();
    int nso       = Process::environment.reference_wavefunction()->nso();
    int nmo       = Process::environment.reference_wavefunction()->nmo();
    int nIrreps   = Process::environment.reference_wavefunction()->nirrep();
    int *orbspi   = Process::environment.reference_wavefunction()->nmopi();
    int *frzcpi   = Process::environment.reference_wavefunction()->frzcpi();
    int *frzvpi   = Process::environment.reference_wavefunction()->frzvpi();
    int *clsdpi   = Process::environment.reference_wavefunction()->doccpi();
    double eNuc   = Process::environment.molecule()->nuclear_repulsion_energy();
    double eOne   = 0.0;
    double eTwo   = 0.0;

    fprintf(outfile, "\n\n\t Irrep  frzcpi doccpi frzvpi  sopi\n");
    for(int h = 0; h < nIrreps; ++h){
        fprintf(outfile, "\t  %3s    %3d    %3d    %3d    %3d\n",
                labels[h], frzcpi[h], clsdpi[h], frzvpi[h], orbspi[h]);
    }
    int nTriMo = nmo * (nmo + 1) / 2;
    int nTriSo = nso * (nso + 1) / 2;
    double *temp = new double[nTriSo];
    Matrix moOei("MO OEI", nIrreps, orbspi, orbspi);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_OEI, temp, nTriMo, 0, 0, outfile);
    moOei.set(temp);
    Matrix moOpdm("MO OPDM", nIrreps, orbspi, orbspi);

    // Before we start, we need to build a mapping array from correlated (no frozen virtuals)
    // to Pitzer, so that the density is ordered consistently with the MO coefficients
    int *toPitzer = new int[nmo];
    size_t corrCount = 0;
    size_t pitzerOffset = 0;
    // Frozen DOCC
    for(int h = 0; h < nIrreps; ++h){
        for(int n = 0; n < frzcpi[h]; ++n){
            toPitzer[corrCount++] = pitzerOffset + n;
        }
        pitzerOffset += orbspi[h];
    }
    // Active DOCC
    pitzerOffset = 0;
    for(int h = 0; h < nIrreps; ++h){
        for(int n = frzcpi[h]; n < clsdpi[h]; ++n){
            toPitzer[corrCount++] = pitzerOffset + n;
        }
        pitzerOffset +=orbspi[h];
    }
    // Active VIRT
    pitzerOffset = 0;
    for(int h = 0; h < nIrreps; ++h){
        for(int n = clsdpi[h]; n < orbspi[h] - frzvpi[h]; ++n){
            toPitzer[corrCount++] = pitzerOffset + n;
        }
        pitzerOffset += orbspi[h];
    }

    if(print > 4)
        for(int n = 0; n < nmo; ++n)
            fprintf(outfile, "Orb %d maps to %d\n", n, toPitzer[n]);

    psio->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
    double **tempOPDM = block_matrix(nmo, nmo);
    psio->read_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *) tempOPDM[0], sizeof(double)*nmo*nmo);
    double *tempMo = new double[nTriMo];
    for(int p = 0; p < nmo; ++p){
      for(int q = 0; q <= p; ++q){
          int P = toPitzer[p];
          int Q = toPitzer[q];
          size_t PQ = INDEX(P,Q);
          tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
      }
    }
    free_block(tempOPDM);
    psio->close(PSIF_MO_OPDM, 1);
    moOpdm.set(tempMo);
    if(print > 4){
        fprintf(outfile, "The MO basis OPDM, in Pitzer order\n");
        print_array(tempMo, nmo, outfile);
        fprintf(outfile, "And again, in matrix form\n");
        moOpdm.print();
        moOei.print();
    }
    delete[] tempMo;

    eOne = eNuc + moOpdm.vector_dot(moOei);

    dpdbuf4 I, G;
    dpd_buf4_init(&G, PSIF_TPDM_PRESORT, 0, ID("[A,A]"), ID("[A,A]"),
                  ID("[A>=A]+"), ID("[A>=A]+"),0, "MO TPDM (AA|AA)");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"),
                  ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    eTwo   = dpd_buf4_dot(&I, &G);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);
    fprintf(outfile, "\n\tMO basis results\n");
    fprintf(outfile, "\tOne electron energy = %16.10f\n", eOne);
    fprintf(outfile, "\tTwo electron energy = %16.10f\n", eTwo);
    fprintf(outfile, "\tTotal energy        = %16.10f\n", eOne + eTwo);

    /*
     * The SO basis
     */
    Matrix soT(nIrreps, orbspi, orbspi);
    Matrix soV(nIrreps, orbspi, orbspi);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_SO_T, temp, nTriSo, 0, 0, outfile);
    soT.set(temp);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_SO_V, temp, nTriSo, 0, 0, outfile);
    soV.set(temp);
    Matrix soOei("SO OEI", nIrreps, orbspi, orbspi);
    soOei.add(soT);
    soOei.add(soV);
    Matrix soOpdm("SO OPDM", nIrreps, orbspi, orbspi);
    psio->open(PSIF_AO_OPDM, PSIO_OPEN_OLD);
    psio->read_entry(PSIF_AO_OPDM, "SO-basis OPDM", (char *) temp, sizeof(double)*nTriSo);
    psio->close(PSIF_AO_OPDM, 1);
    soOpdm.set(temp);

    if(print > 4){
        soOpdm.print();
        soOei.print();
    }
    eOne  = eNuc + soOpdm.vector_dot(soOei);

    dpd_buf4_init(&G, PSIF_TPDM_HALFTRANS, 0, ID("[n,n]"), ID("[n,n]"),
                  ID("[n>=n]+"), ID("[n>=n]+"), 0, "SO Basis TPDM (nn|nn)");
    dpd_buf4_init(&I, PSIF_SO_PRESORT, 0, ID("[n,n]"), ID("[n,n]"),
                  ID("[n>=n]+"), ID("[n>=n]+"), 0, "SO Ints (nn|nn)");
    eTwo   = dpd_buf4_dot(&I, &G);
    dpd_buf4_close(&I);
    dpd_buf4_close(&G);

    psio->close(PSIF_TPDM_HALFTRANS, 1);
    psio->close(PSIF_LIBTRANS_DPD, 1);
    psio->close(PSIF_TPDM_PRESORT, 1);
    psio->close(PSIF_SO_PRESORT, 1);
    fprintf(outfile, "\n\tSO basis results\n");
    fprintf(outfile, "\tOne electron energy = %16.10f\n", eOne);
    fprintf(outfile, "\tTwo electron energy = %16.10f\n", eTwo);
    fprintf(outfile, "\tTotal energy        = %16.10f\n", eOne + eTwo);
    fflush(outfile);
    delete [] temp;

    return Success;   
}

}} // End Namespace