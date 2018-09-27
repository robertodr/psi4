import psi4

psi4.set_output_file('output.dat', False)

h2o = psi4.geometry("""
  0 1
  H
  O 1 0.957
  H 2 0.957 1 104.5
""")

psi4.set_options({
  'freeze_core': 'false',
  'basis': 'sto-3g',
  'scf_type': 'pk'
})

scf_e, scf_wfn = psi4.energy('SCF', return_wfn=True)
psi4.core.print_variables()
psi4.compare_values(-74.9629054164371809, scf_e, 6, 'SCF PK Energy')

cc_wfn = psi4.core.CCWavefunction(scf_wfn)
ccsd_e = cc_wfn.compute_energy()
psi4.core.print_out('CCSD energy {:20.12f}'.format(ccsd_e))

# This is what we want to emulate
#ccsd_e, wfn = psi4.properties('ccsd',properties=['dipole'],return_wfn=True)
#psi4.oeprop(wfn,"DIPOLE", "QUADRUPOLE", title="(OEPROP)CC")
#psi4.compare_values(psi4.get_variable("(OEPROP)CC DIPOLE X"), 0.000000000000,6,"CC DIPOLE X")             #TEST
#psi4.compare_values(psi4.get_variable("(OEPROP)CC DIPOLE Y"), 0.000000000000,6,"CC DIPOLE Y")             #TEST
#psi4.compare_values(psi4.get_variable("(OEPROP)CC DIPOLE Z"),-1.840334899884,6,"CC DIPOLE Z")             #TEST
#psi4.compare_values(psi4.get_variable("(OEPROP)CC QUADRUPOLE XX"),-7.864006962064,6,"CC QUADRUPOLE XX")   #TEST
#psi4.compare_values(psi4.get_variable("(OEPROP)CC QUADRUPOLE XY"), 0.000000000000,6,"CC QUADRUPOLE XY")   #TEST
#psi4.compare_values(psi4.get_variable("(OEPROP)CC QUADRUPOLE XZ"), 0.000000000000,6,"CC QUADRUPOLE XZ")   #TEST
#psi4.compare_values(psi4.get_variable("(OEPROP)CC QUADRUPOLE YY"),-4.537386915305,6,"CC QUADRUPOLE YY")   #TEST
#psi4.compare_values(psi4.get_variable("(OEPROP)CC QUADRUPOLE YZ"), 0.000000000000,6,"CC QUADRUPOLE YZ")   #TEST
#psi4.compare_values(psi4.get_variable("(OEPROP)CC QUADRUPOLE ZZ"),-6.325836255265,6,"CC QUADRUPOLE ZZ")   #TEST

psi4.core.print_variables()
