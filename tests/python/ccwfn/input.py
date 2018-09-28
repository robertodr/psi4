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
  'scf_type': 'pk',
  'wfn' : 'ccsd'
})

scf_e, scf_wfn = psi4.energy('SCF', return_wfn=True)
psi4.core.print_variables()
psi4.compare_values(-74.9629054164371809, scf_e, 6, 'SCF PK Energy')

psi4.core.print_out('CCWavefunction\n')
cc_wfn = psi4.core.CCWavefunction(scf_wfn)
ccsd_e = cc_wfn.compute_energy()
psi4.core.print_out('CCSD energy {:20.12f}\n'.format(ccsd_e))
psi4.core.print_variables()

# This is what we want to emulate
ccsd_e, wfn = psi4.energy('ccsd', return_wfn=True)
psi4.core.print_variables()
