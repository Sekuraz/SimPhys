     &SCF
         ADDED_MOS 50
         &DIAGONALIZATION
            ALGORITHM STANDARD
            EPS_ADAPT 0.01
         &END DIAGONALIZATION
         &SMEAR  ON
            METHOD FERMI_DIRAC
            ELECTRONIC_TEMPERATURE [K] 300
         &END SMEAR
         &MIXING
            METHOD BROYDEN_MIXING
            ALPHA 0.2
            BETA 1.5
            NBROYDEN 8
         &END MIXING
      &END SCF
      &KPOINTS
         SCHEME MONKHORST-PACK 2 2 2
         WAVEFUNCTIONS REAL
         SYMMETRY .FALSE.
         FULL_GRID .FALSE.
         PARALLEL_GROUP_SIZE -1
      &END KPOINTS
      &PRINT
         &PDOS
            FILENAME ./diamond_pdos
            NLUMO -1
         &END PDOS
      &END PRINT


