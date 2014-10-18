Attached are the source and pipeline for installing the latest version of the PaCE software.

(i) checksum - for verifying the file downloads (with md5sum),
(ii) PaCE_release - has the source code which needs to be compiled using the make command (see PaCE-Build-Readme),
(iii) PaCE-pipeline.2.2 - has the directory structure from which to run the PaCE program

SYSTEM REQUIREMENTS

Any linux/unix based cluster with MPICH installed on it.

INSTALLATION STEPS

1. First verify checksum
2. Build the PaCE binary from the PaCE_release folder using the make command. Depending on the mpich location on the cluster environment, the Makefile has to be changed before running the make command.
3. Copy the "PaCE" binary to PaCE-pipeline2.2/PaCE-clustering/
4. Copy "PreprocessPaCE.pl"  to PaCE-pipeline2.2/PaCE-clustering/

TESTING

From PaCE-clustering folder, run the following commands:

- PreprocessPaCE.pl ../datafiles/tEST.data 0 > ../datafiles/tEST.data.PaCE
- mpirun -np 4 PaCE ../datafiles/tEST.data.PaCE 500 PaCE.cfg

If there were no errors at any stage then it means that the installation was a success. 
