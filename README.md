# Exact Ligand Solid Angles

## INPUT
**ComplexDataBase1.txt**

Cartesian coordinates for the complex are stored in ComplexDataBase1.txt in the following format. *The apex atom (central atom) must be listed first!*
```
(XN) NAME
1 Z x y z
2 Z x y z
0
```
where N is the index of the complex in the database (to be called upon by SolidAngleDriver.nb), NAME is a unique identifier (for you to identify, it's not used in the program), Z is the atomic symbol (see SolidAnglePackage.nb section for list of supported atoms), and x, y, and z are the Cartesian coordiantes of the corresonding atom in Å. 0 follows each input to separate complexes.


**ComplexDataBase.nb/ComplexDataBase.m**
ComplexDataBase.nb reads in the Cartesian coordinates and formats Apex, Ligands, and XAtoms correctly for SolidAngleDriver.nb. To change the database .txt file where the complexes are stored, edit the *lingandsIn* variable. For example, for a new database NewLigands.txt, the *lingandsIn* variable should be ligandsIn=OpenRead[“NewLigands.txt”]. To solve the solid angles of multiple ligands on a single apex, edit the variable *Ligands*. For example, if atom 1 is the metal center, atoms 2-6 belong to the first ligand, and atoms 7-14 belong to the second ligand, the variable *Ligands* sould be Ligands={Range[2,6],Range[7,14]}. Once edits have been made to ComplexDataBase.nb, the file must be saved as a Mathematica Package with the .m extension to allow SolidAngleDriver.nb to call upon the package.


## RUN FILES
**SolidAngleDriver.nb**
SolidAngleDriver.nb is the only package the user needs to modify to use the code as is. This interface calls on ComplexDataBase.m to read Cartesian coordinates from ComplexDataBase1.txt and calls SolidAnglePackage.m to calculate the solid angle and visualize the results. The variable *SolidAngleDirectory* should be defined as the directory with the ExactLigandSolidAngle files. *ComplexSet* defines the complex(es) of interest. The index should correspond to N in ComplexDataBase1.txt. For example, to find the solid angles of complexes 12-17, ComplexSet = Range[12, 17] would be used. To find the solid angle of only complex 12, the range should be ComplexSet=Range[12,12].

The output can be controlled by the variable *kPrint*:
* *kPrint* = 0, No printing within package
* *kPrint* = 1, Print total solid angle and dissection of loop contributions for each ligand
* *kPrint* = 2, For each ligand print a (non-shaded) 3D plot highlighting the border arcs and showing the shadow cones intersecting with the unit sphere
* *kPrint* = 3, For each ligand print a table of atomic radii, vertex angles, and Cartesian coordinates
* *kPrint* = 4, For each ligand print the aforementioned 3D plot with labels on the arcs and for comparison also tabulate the dissection of loop contributions for each ligand
* *kPrint* = 5, For each ligand print a 3D plot showing the border arcs and filling gray into the shadow cones intersecting with the unit sphere

*Alternate Atomic Radii.* The choice of atomic radii can be changed with the variable *kR*. When kR = 1 the atomic radii from [Bondi, A. J. Phys. Chem. 1964, 68, 441-451.] will be used. When kR = 2 the zero energy point radii from [Guzei, I. A.; Wendt, M. Dalton Trans. 2006, 3991-3999] will be used. As of now only H, He, C, N, O, F, Ne, Si, P, S, Cl, Ar, As, Se, Br, Kr, Te, I, and Xe are included. Radii for additional atoms can be added in the SolidAnglePackage.m file. The implemented Bondi atomic radii are {"H", "He", "C", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe"} with corresponding values {1.20, 1.40, 1.70, 1.55, 1.52, 1.47, 1.54, 2.10, 1.80, 1.80, 1.75, 1.88, 1.85, 1.90, 1.85, 2.02, 2.06, 1.98, 2.16} in Å. The implemented zero energy point radii are {"H", "He", "C", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe"}, with corresponding values {1.000, 1.311, 1.539, 1.521, 1.470, 1.413, 1.350, 1.834, 1.801, 1.757, 1.599, 1.649, 1.879, 1.861, 1.845, 1.831, 1.955, 1.941, 1.928} in Å.

*Closing Gaps.* If a ligand profile contains a small indentation, that gap may not be relevant in quantifying the total sterics. We added an option to effectively close that gap by placing a dummy atom that overlaps with the two atoms at the start of the gap. The option is controlled by setting the variable *PlugGap* to True or False. To determine the two atoms to use, first find the solid angle with *kPrint* set to 4 and *PlugGap* set to False. The 3D plot will have each arc labeled. From the labeling, determine which two atoms should be plugged to fill the gap. To find the solid angle with a plugged gap, set PlugGap = True, arc1 = number of first arc in plug, arc2 = second arc in plug, and plugradius = radius of the plug. The value of *plugradius* must be determined manually by checking the 3D profile to confirm that the gap is closed.


**SolidAnglePackage.nb/SolidAnglePackage.m**
SolidAnglePackage.nb includes all mathematics used to solve for and visualize the solid angle. A list of standard atomic radii and zero energy point radii are already written into the program but can be changed or added to. To modify the list of atomic radii, locate the variable *RvdW* (Bondi atomic radii) and *RZ* (zero energy point radii). *RvdW* are *RZ* are lists of three lists. The first list holds atomic symbols whose position in this list correspond with values given in the next two lists, the second list holds atomic radii in Å, and the third list holds the atom colors for plotting. Once edits have been made to SolidAnglePackage.nb, the file must be saved as a Mathematica Package with the .m extension to allow SolidAngleDriver.nb to call upon the package.


## OUTPUT
Program output will appear in SolidAngleDriver.nb along with the total run time. Output level is controlled by *kPrint* in SolidAngleDriver.nb.


## CITATION
If you find this work useful, please cite our publication: 
* [Jenna A. Bilbrey  Arianna H. Kazez  Jason Locklin  Wesley D. Allen, "Exact Ligand Solid Angles", Journal of Chemical Theory and Computation, 2013, 9 (12), pp. 5734–5744.](https://pubs.acs.org/doi/abs/10.1021/ct400426e)

