
Steps to obtain the nightly tests for the band-to-trap tunneling model: 

1. Run npn-hbt.jou using Cubit to obtain npn-hbt.exo (Depending on the Cubit version, there may be an error saying some of the IDs are not found. Solution: run the journal file step by step, find the new IDs, and replace the old IDs with the new ones.)

2. Run "epu -decomp 4 npn-hbt.exo" to obtain the parallel exodus input file.

3. Run "mpirun -np 4 charon_mp.exe --i=npn-hbt.nlp.xml" to obtain the nlp solution. 

4. Run "mpirun -np 4 charon_mp.exe --i=npn-hbt.dd.0V.telt.xml --current" to obtain the equilibrium solution at 0 V. 

5. Run "mpirun -np 4 charon_mp.exe --i=npn-hbt.dd.telt.xml --current" to obtain a gummel sweep (Vbe = 0 to 1.35 V) without the Trap SRH model. The current data are saved in current-dd-telt.dat.

Can also run "mpirun -np 4 charon_mp.exe --i=npn-hbt.sg.telt.xml --current" for the SGCVFEM discretization. The current data are saved in current-sg-telt.dat. 

6. Run "mpirun -np 4 charon_mp.exe --i=npn-hbt.dd.telt.0.8V.xml --current" to get the solution at Vbe = 0.8 V without the Trap SRH model. 

Can also run "mpirun -np 4 charon_mp.exe --i=npn-hbt.sg.telt.0.8V.xml --current" for the SGCVFEM discretization.

7. In the subfolders: dd-schenk, dd-newdos, sg-schenk, and sg-newdos, the "Schenk ConstFDOS" and "Schenk NewDOS" tunneling models are tested on a 1D InGaP/GaAs HBT (from the literature) for the FEM-SUPG and CVFEM-SG discretization methods. 

The "Schenk ConstFDOS" tests are set up to run for 6 voltage points, while the "Schenk NewDOS" tests only run for 2 voltage point to save run time for nightly tests. The "Schenk NewDOS" model is also tested under the heavytest directory for a LOCAL sweep from Vbe=0.8 V to 1.3 V.   


Other non-nightly tests included in the folder: 

8. Run "mpirun -np 4 charon_mp.exe --i=npn-hbt.dd.trap.nob2t.xm --current" to obtain a gummel sweep (Vbe = 0.8 to 1.3 V) with the Trap SRH model but without the band-to-trap tunneling. The simulation uses the solution at Vbe=0.8 V from Step 5 as the initial guess. The current data are saved in current-dd-trap-nob2t.dat. 

Can also run "mpirun -np 4 charon_mp.exe --i=npn-hbt.sg.trap.nob2t.xml --current" for the SGCVFEM discretization. The current data are saved in current-sg-trap-nob2t.dat.


