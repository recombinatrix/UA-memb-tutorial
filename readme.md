# Tutorial: Gromacs United Atom System Setup

This was written by Ada Quinn, a researcher in computational chemisty at The University of Queensland.  This readme is version 1.1.2, dated 18 March 2023.  It is distributed under a GNU GPL-v3 (2007) licence.

I wrote this to provide to students in my own research group, where I would be on hand to support them if something went wrong. If you are not my student please know that I hope it is helpful to you, but use it at your own risk.

## Introduction

This is a tutorial for building a simulation system to study the human glycine transporter GlyT2 (uniprot id [SLC6A5](https://www.uniprot.org/uniprotkb/Q9Y345/entry)).  It was written to work with Gromacs version 2021.4, using the [GROMOS 54a7](https://doi.org/10.1007/s00249-011-0700-9) united atom forcefield.

In this tutorial, you will place a homology model of GlyT2 into a binary POPC/cholesterol membrane, solvate this protein/membrane system, introduce sodium and chloride ions at a physiologically releveant concentration, and perform equilibration.  

This entire tutorial and all associated files is availible for download from github at [https://github.com/recombinatrix/UA-memb-tutorial](https://github.com/recombinatrix/UA-memb-tutorial).  If you have git installed on your computer, you can download the tutorial and all associated files by running

~~~s
git clone https://github.com/recombinatrix/UA-memb-tutorial.git
~~~

### Requirements

This tutorial assumes you are working with linux and have basic linux knowledge.

You will need working installation of gromacs 2021.4, visualising molecular dynamics (VMD), and PyMOL.  You are assumed to know the basics of using VMD.  If you need to learn how to use VMD, I recommend working through [these tutorials](https://www.ks.uiuc.edu/Training/Tutorials/vmd-index.html) from the [Theoretical and Computational Biophysics Group at the University of Illinois Urbana-Champaign](https://www.ks.uiuc.edu/)

### Useful resources

The [Gromacs manual](https://manual.gromacs.org/2021.4/index.html) includes extensive documentation for every gromacs command used in this tutorial.  Whenever you use a new gromacs command, you should have a look at the manual to see what it does, and see if you can figure out why I wrote it in that particular way.

I strongly encourage you to learn to work with an integrated development environment, such as [Visual Studio Code](https://code.visualstudio.com/).  IDEs allow you to explore your file system, edit multiple documents at once, and run multiple terminal sessions all in a single window.  There are helpful extensions that make allow vscode to understand the shape of gromacs file types, which makes them significantly easier to read and understand.  This is especially valuable while you are learning how to work with gromacs and the associated file formats for the first time.

I suggest that you install VSCode and a [gromacs helper extension](https://marketplace.visualstudio.com/items?itemName=SupernovaZJX.gmx-helper) before you start this tutorial.

### Recordkeeping

Whenever you use a computer to perform scientific process, such as building or analysing a molecular dynamics simulation, you ***must*** make a readme file.  A readme is a written record of everything you did, including every command you ran and every file you changed, along with a minimal explanation as to what you were trying to accomplish.  

Your readme file should be sufficiently detailed that another scientists (who does not know you and does not know what you intended to do) could sit down at your computer and replicate your work with only your readme file.  Make sure your readme is good, because sometimes that scientist who does not know you, and does not know what you intended to do, will be you in the future, trying to figure out what you did six months ago!  

You need to make your own readme file while you work through this tutortial.  **Do not try to write your readme file after you have finished.  That way lies misery.**

## The Tutorial

### Getting started

Any simulation needs three things:

* A set of *instructions*, which tells gromacs how to perform the simulation.  These are usually contained in `.mdp` and `.sh` files, and inside the gromacs code.
* A set of *coordinates*, which describes the position of every atom in the system.  For the kind of work I do, these are usually contained in a `.pdb` or `.gro` file.
* A *forcefield* and associated *parameters*, which describes the forces between atoms, and how different atoms are connected.   These are usually contained in a `.top` file, a number of `.itp` files, and a large number of files contained in a `.ff` folder.

For your tutorial, some of these things are provided.

The *instructions* for running each of the molecular dynamics steps have been provided through a series of `.mdp` and `.sh` files.

The *coordinates* for GlyT2 and a POPC/CHOL membrane have been provided for you.  Over the course of this tutorial you will need to combine them into a single set of coordinates.

`GlyT2_capped.pdb` contains coordinates of a homology model of GlyT2 based on the homologue dDAT, [published here](https://doi.org/10.1371/journal.pone.0157583).  `POPC_CLR_550.gro` contains coordinates for a model membrane of 80% POPC and 20% Cholesterol, made with [memgen](https://memgen.uni-saarland.de/), and equilibrated by me.

The *forcefield* you will use is the [GROMOS 54a7](https://doi.org/10.1007/s00249-011-0700-9) forcefield, which has been provided in the `gromos54a7.ff/` folder, along with some useful `.itp` files. The forcefield folder I have provided includes instructions for the forcefield physics, and descriptions of some specific molecules like lipids and amino acids (*parameters*).  We will also need to generate *parameters* for additional molecules we introduce into our system, like our GlyT2 protein.

### Generating GlyT2 parameters

Lets start by generating parameters for GlyT2.

We generate a protein topology using `pdb2gmx`.

~~~s
gmx pdb2gmx -f GlyT2_capped.pdb -vsite hydrogens -heavyh -ter -ignh -o 01_GlyT2_capped_pdb2gmx.pdb -p GlyT2_POPC_CLR.top
~~~

Choose these options:

* The `GROMOS54a7` forcefield folder in the folder (probably option 1)
* The `SPC` water model (option 1)
* No change to the N-terminal cap: `None` (option 2)
* No change to the C-terminal cap: `None` (option 2)

This will generate two files: a coordinate file called ` 01_GlyT2_capped_pdb2gmx.pdb ` , which contains all atoms in the correct positions, and a topology file ` GlyT2_POPC_CLR.top `, which describes the parameters of GlyT2, and calls the GROMOS 54a7 forcefield.

First, have a look at ` GlyT2_POPC_CLR.top `.  At the beginning is a series of comments that describe how , where, and when it was made.

Next there are some rules for the forcefield, like

~~~c
    #define HEAVY_H
~~~

which tells gromacs to use the rules for heavy hydrogens, and

~~~c
    #include "gromos54a7.ff/forcefield.itp"
~~~

which tell gromacs to use the GROMOS 54a7 force field from the file ` gromos54a7.ff/forcefield.itp `.

Next, there is a *very* long section that described the topology of the protein.  This includes the charge, mass and atomtype of every atom, as well as every defined bond, every defined bond angle, and every defined dihedral, as well as various other interactions and constraints.

This section starts with:

~~~c
    [ moleculetype ]
    ; Name            nrexcl
    Protein             3
~~~

and ends with

~~~c
    5254  5248  5249  5250     4 
    5397  5392  5393  5394    -4 
    5398  5392  5393  5394     4 

    ; Include Position restraint file
    #ifdef POSRES
    #include "posre.itp"
    #endif
~~~

After that, we have some rules to include water and ion topology, again calling `.itp` files for water and ions.

Right at the end, there is a section that lists every molecule in the system.  Right now, that's very simple because there's only one molecule: the protein.

~~~c
    [ molecules ]
    ; Compound        #mols
    Protein             1
~~~

The topology file is very big and unwieldy, you will need to edit it.  Usually we like to keep topology files quite minimal, and use ` .itp ` include files to contain the properties of individual molecules within the system.  This is very useful when you need to build your next system, as you can just choose the `.itp` files you need.

Now, make a new file called `GlyT2.itp`.  Copy all of the GlyT2 topology
information from your `.top` file (`GlyT2_POPC_CLR.top`) into the `.itp` file.  You should be copying everything from

~~~c
    [ moleculetype ]
    ; Name            nrexcl
    Protein             3
~~~

to

~~~c
    ; Include Position restraint file
    #ifdef POSRES
    #include "posre.itp"
    #endif
~~~

Then, delete that entire section of the `.top` file, and replace it with this line:

~~~c
    #include "GlyT2.itp"
~~~

You've just put all that information into a single include file

### Energy minimization

We need to define the size of our simulation system, by specifying the box size of our universe.  We will do this with `editconf`.

Before you run this command, go to the [gromacs manual](manual.gromacs.org) and look at the entry for `editconf`.  Make sure you are looking at the manual for the version of gromacs you are using!

What does this `editconf` do?  What do the `-f`, `-o` and `-box` flags do?

~~~s
gmx editconf -f 01_GlyT2_capped_pdb2gmx.pdb -o 02_GlyT2_capped_pdb2gmx_box.pdb \
                 -box 17 17 12
~~~

Next, do an energy minimization on this system.  This is a short simualtion that tries to move the atoms in the simualtion towards the nearest local minima in the energy landscape.  This is great for removing artefacts from the way you assemble your system, such as atoms with overlapping van der waals radii.  First we will use `grompp`, the gromacs preprocessor, to assemble all the parts of your simulation into a single tpr file, and perform some quality checks to make sure the simulation is well formed. Then we will use `mdrun` to actually run the energy minimzation.

Look at the gromacs manual for `grompp`.  What does `grompp` stand for?  What does `grompp` do?  Which flags are you using, and what do they do?

~~~s
gmx grompp -f minimise.mdp -c 02_GlyT2_capped_pdb2gmx_box.pdb                  \
           -p GlyT2_POPC_CLR.top -o 03_GlyT2_capped_EM.tpr -maxwarn 1
~~~

Try running grommp without the `-maxwarn 1` flag, and see what happens.

~~~s
gmx mdrun -v -deffnm 03_GlyT2_capped_EM
~~~

Look at the output and see if the numbers are sensible. Every energy minimization step in this tutorial should end with a negative value for potential energy, and with the largest force on any one atom smaller than 10^4

### Put the protein in bilayer

GlyT2 is a membrane protein, so we need to add a bilayer.  There is a POPC/CLR bilayer in the file ` POPC_CLR_550.gro `

This file has water and ions in them, and we want to get rid of those for now.  Open ` POPC_CLR_550.gro ` in  VMD, and save the coordinates of just the lipid bilayers, using the selection string ` resname POPC CLR `.  Save this file as ` 04_POPC_CLR_550_dry.pdb `  Open this file in vmd to make sure it looks correct.

Now you need to position the protein so it is overlapping the membrane in the correct orientation.

04_POPC_CLR_550_dry.pdb is very short in the Z dimension, and there isn't enough room for GlyT2.  We need to put the membrane in a bigger box so we can arrange the protein and membrane in the correct orientation.  To make a bigger box, use `editconf`

~~~s
gmx editconf -f 04_POPC_CLR_550_dry.pdb -o 05_POPC_CLR_550_dry_box.pdb -box 17 17 12
~~~

Now we need to position the protein in the membrane, in the correct orientation.  To do this, we need to know what the correct orientation is!  We're going to use the orientation of a related protein, dDAT, as a reference.  You can download to orientation dDAT in a membrane from [the OPM database entry 4m48](https://opm.phar.umich.edu/proteins/2264), by clicking ` Download OPM File: 4m48.pdb `

Open ` 05_POPC_CLR_550_dry_box.pdb ` and ` 4m48.pdb ` in VMD.  The OPM structure has a red and blue layers that correspond to the top and bottom of a membrane.  Try to line these up with the top and bottom of ` 05_POPC_CLR_550_dry_box.pdb `, as well as you can.  You can use the phosphate atoms from POPC headgroups as a proxy for the interface of the bilayer, with the selection ` name P ` and the draw method `VDW`.

To move 4m48 up and down you can either use `gmx editconf` or the vmd movement tools.  See what works best for you.  

Once you have the protein at the correct height we need to mvoe it to the center of the membrane patch. Look from the top, and move the protein around until it sits in the center of the membrane.  Look from the side again to check it's still at the correct height.  Once you are satisfied, save it as ` 05_4m48_aligned.pdb `.

Next we want to align our GlyT2 model with 4m48.  To do this we're going to use PyMOL.  

Open`05_4m48_aligned.pdb` and `03_GlyT2_capped_EM.gro` in PyMOL.  In the right hand menu bar, click on the `A` next to `03_GlyT2_capped_EM.gro`, and choose `align > to molecule > 05_4m48_aligned.pdb`   **NB: Make sure you align GlyT2 to 4m48, not the other way around!**

The two proteins should now be very closely aligned, with their helices almost overlapping. Now we need to save th new coordinates.  Go to `File > Export Molecule`.  For selection choose ` 03_GlyT2_capped_EM `.  Check the box that says ` original atom order (according to "rank") ` (so PyMOL doesn't rearrange the atoms) and click `save…`. In the save dialogue box choose the `pdb` file format (NOT the default, which is usually `PDBx/mmCIF` ).  Save the data as ` 05_GlyT2_capped_moved.pdb `

Now combine the aligned protein and the dry membrane into one file with `cat`

~~~s
cat 05_GlyT2_capped_moved.pdb 05_POPC_CLR_550_dry_box.pdb > 06_GlyT2_POPC_CLR_crude.pdb
~~~

Open ` 06_GlyT2_POPC_CLR_crude.pdb `.  Search for the place in the middle where the two files are joined.  It'll start with a line like ` TER ` or ` END `, and it'll look something like this:

~~~pdb
    …
    ATOM   5468  O   ASP   757      26.817   8.144 -29.045  1.00  0.00           O
    TER
    ATOM   5469  N   NH2   758      24.674   7.205 -29.068  1.00  0.00           N
    ATOM   5470  H1  NH2   758      24.122   6.371 -29.061  1.00  0.00            
    ATOM   5471  H2  NH2   758      24.228   8.100 -29.097  1.00  0.00            
    TER
    ENDMDL
    CRYST1  169.833  169.833   84.540  90.00  90.00  90.00 P 1           1
    ATOM      1  CN1 POPCX   1       4.670 162.650  64.350  0.00  0.00            
    ATOM      2  CN2 POPCX   1       6.150 164.210  63.380  0.00  0.00            
    ATOM      3  CN3 POPCX   1       4.110 164.890  64.570  0.00  0.00            
    …
~~~

get rid of every line that doesn't start with ` ATOM ` from this middle part, so it looks something like this:

~~~pdb
    …
    ATOM   5468  O   ASP   757      26.817   8.144 -29.045  1.00  0.00           O
    ATOM   5469  N   NH2   758      24.674   7.205 -29.068  1.00  0.00           N
    ATOM   5470  H1  NH2   758      24.122   6.371 -29.061  1.00  0.00            
    ATOM   5471  H2  NH2   758      24.228   8.100 -29.097  1.00  0.00            
    ATOM      1  CN1 POPCX   1       4.670 162.650  64.350  0.00  0.00            
    ATOM      2  CN2 POPCX   1       6.150 164.210  63.380  0.00  0.00            
    ATOM      3  CN3 POPCX   1       4.110 164.890  64.570  0.00  0.00            
    …
~~~

**NB:  Only edit the middle!**  Don't delete any lines from the very beginning or very end of the file!

Save your edits and close the file.

#### Identifying problems

Open ` 06_GlyT2_POPC_CLR_crude.pdb ` in VMD.  There are two problems with this file:

To find the first problem, make a representation for the POPC lipids, the CLR lipids, and the protein.  For the protein, you need to use a special selection string:  `protein or resname ACE NH2`.  This is because sometimes VMD does not recognise the `ACE` and `NH2` caps as part of the protein.

Look at your representations.  You can see that the protein and the membrane are overlapping; some of the lipids are in the same place as GlyT2. We need to fix this by removing every lipid that is overlapping with the protein.

The second problem is more subtle.  To see it, we're going to need our representation for to show the protein: `protein or resname ACE NH2`.  Next, we need a second representation:  `same resid as (protein or resname ACE NH2)`.  What do you see?  Try making a third representation:  `(same resid as (protein or resname ACE NH2)) and not (protein or resname ACE NH2)`.  Compare it to the previous two.  

What do you see?  Stop and think about what these representations are doing, and why they work differently.  Why is `same resid as (protein or resname ACE NH2)` not the same as `(same resid as (protein or resname ACE NH2)) and not (protein or resname ACE NH2)`?

What is the second problem?  We need to fix this problem first.

#### Fixing the problem with the membrane residue numbers

To fix our membrane, we need to renumber the lipids.  Go to the [gromacs manual](manual.gromacs.org), and find the entry for the `editconf` command.  You need to write a command to change membrane residue numbers.  Look at the options in the manual and find a way to change the residue numebrs.  Think about what the residue number of the first POPC lipid needs to be.  

Run your command to renumber the lipids.  Call your new membrane `06_POPC_CLR_550_dry_renumbered.pdb`.  Look at the file and see if it worked correctly.

Combine your renumbered membrane with GlyT2, and clean it up by removing the lines between the protein and the membrane:

~~~s
cat 05_GlyT2_capped_moved.pdb 06_POPC_CLR_550_dry_renumbered.pdb > 06_GlyT2_POPC_CLR_crude_renumbered.pdb
~~~

#### Fixing the problem with the overlapping lipids

Open `06_GlyT2_POPC_CLR_crude_renumbered.pdb` in vmd.

We need to remove the lipids overlapping with GlyT2.  We're going to do this by finding a selection string for all the atoms we want to keep, and using it to save the coordinates of these atoms.

To start with, lets think about some selection strings.

1. `(resname POPC CLR and within 5 of protein)`.  This gives you every atom that meets two conditions: One, the atom belongs to a POPC or CLR molecule.  Two: the atom is within 5 A of one of the atoms in GlyT2.  

2. `not (resname POPC CLR and within 5 of protein)`.  This gives you every atom that doe not meet the two conditions above.  It's the opposite of the previous selection.

3. `(not (protein or resname NH2 ACE) and not (resname POPC CLR and within 5 of protein)`.  This is the same as above, but we have added a new condition:  The atoms must also not be part of the protein.

How is number three different to number two?  Why have we added the extra condition?

Number three is still not good enough.  Look at the POPC molecules near GlyT2.  You can see some of them have been cut in half. This is no good.

Think about the VMD selections you have learned in this section, and in the  **Identifying problems** section, and use VMD to create a file that contains GlyT2, and and the lipids you want to keep, and none of the lipids you want to remove.  Make sure you don't cut any lipids in half.  Make sure you don't remove any parts of the protein by accident.

Some people like to do this with just one selection string.  Some people like to do this by saving the protein to one file, saving the lipids they want to keep to another file, and those files afterwards combining them together. Find a process that works for you.  Save your final file, with GlyT2 and the lipids, as ` 07_GlyT2_POPC_CLR_hole.pdb `  If you do it wrong, you'll know because your next `grompp` will fail.

Open `07_GlyT2_POPC_CLR_hole.pdb` in vmd and check it looks correct.  There should be a protein in a membrane, and there should be a big hole in the membrane around the protein.  There should be no lipids touching or overlapping with the protein.  If this is not what you see, keep trying until you get it right.  `hole.tga` contains an example of what it should look like from the top.

Once you are happy you have everything right, you need to edit your topology file to include the new lipids.

Edit ` GlyT2_POPC_CLR.top ` to add `#include` statements for the `.itp` files for popc and cholesterol.  Have a look at ` POPC_CLR_550.top ` to see what it should look like.

Next, count how many POPC molecules are present.  To do this, you can use `grep`.  Grep looks for particular patterns in a file.

If we use the command

~~~s
grep -c " CA " 07_GlyT2_POPC_CLR_hole.pdb
~~~

it will tell us how many times the pattern "`CA`" (with spaces) occurs in your file.  We could use this to count alpha carbons.  Unfortunately, alpha carbons are found in both POPC and amino acids, so that won't do.  If you ran

~~~s
grep -c " POPC " 07_GlyT2_POPC_CLR_hole.pdb
~~~

it would count the number of times "`POPC`" occurs in the file (again, with spaces).  Unfortunately, the string "`POPC`" occurs multiple times in every POPC molecule, so that won't do either.

Open `07_GlyT2_POPC_CLR_hole.pdb` as in a text editor, and scroll down to find a single POPC molecule.  See if you can come up with some string of characters that occurs exactly once in every single POPC molecule, and nowhere else.

Add the count of POPC molecules to the end of `GlyT2_POPC_CLR.top`

It should look something like this

~~~c
    [ molecules ]
    ; Compound        #mols
    Protein             1
    POPC               880
~~~

Next, come up with a pattern to count the cholesterol molecules, and add them as well. You'll need to figure this one out yourself. Once again, you want a string of characters that occurs exactly once in every cholesterol moelcule, but nowhere else.

Save the changes you made to your `.top` file, and do an energy minimisation.  If your `grompp` fails, you've probably made a mistake in one of the previous steps.  Have a look at your error message and see if you can figure out where.

~~~s
gmx grompp -f minimise.mdp -c 07_GlyT2_POPC_CLR_hole.pdb -p GlyT2_POPC_CLR.top \
           -o 08_GlyT2_POPC_CLR_hole_EM.tpr -maxwarn 1

gmx mdrun -v -deffnm 08_GlyT2_POPC_CLR_hole_EM
~~~

Hopefully this worked.    Remember to check your evergy levels.

### Solvation

Next, we're going to add solvent to our system.

We can do this with the ` gmx solvate ` command.  But first, look up `solvate`  in the gromacs manual.  What does it do?  What flags are you using?

~~~s
gmx solvate -cp 08_GlyT2_POPC_CLR_hole_EM.gro -p GlyT2_POPC_CLR.top           \
            -o 09_GlyT2_POPC_CLR_water_crude.pdb
~~~

gromacs will add solvent to your system, and update your `.top` file automatically.  whenevr gromacs overwrites a file, it saves a backup, and there will be a backup of your original `.top` file called something like `#GlyT2_POPC_CLR.top.1#`.  **IMPORTANT:** other programs, like `cat` or `cp` or `mv` will **NOT** make a backup of your files.  Often, they will just overwrite your files immediately and forever, without even asking.  Always be careful before you run a command that will alter or move a file.

Open ` 09_GlyT2_POPC_CLR_water_crude.pdb ` in vmd.  You can see there is a lot of water that has ended up inside the bilayer.  This is because `gmx solvate` doesn't know we're using a united atom system, and it thinks there's a lot of empty space inside our membrane. You need to remove this errant solvent manually.

remember you can visualise the interface of the bilayer with the selection ` name P `, which shows the phosphate atoms at the top and bottom of the membrane.

Use vmd selections to find the water that is inside the bilayer.  It will be something like ` resname SOL and (same residue as resname SOL and z < 100 and z > 50) `.  Change the numbers until you find the correct heights so no water is between phosphate atoms.  

Save the coordinates of everything except that water.  You'll want a selection string like ` not ( resname SOL and (same residue as resname SOL and z < 100 and z > 50) ) `.  Save this file as ` 09_GlyT2_POPC_CLR_water_cut.pdb `

Open it in VMD and make sure it looks correct.  Once you're happy with it, figure out how you can use `grep -c` to count the number of remaining water molecules.  Update your topology file to reflect the new number.  Then, do another energy minimzation.

~~~s
gmx grompp -f minimise.mdp -c 09_GlyT2_POPC_CLR_water_cut.pdb -p GlyT2_POPC_CLR.top -o 10_GlyT2_POPC_CLR_water_EM.tpr -maxwarn 1

gmx mdrun -v -deffnm 10_GlyT2_POPC_CLR_water_EM
~~~

If your `grompp` fails, you've probably made a mistake in one of the previous steps.  Have a look at your error message and see if you can figure out where.

GlyT2 exists in cell membranes, and in a real physiological system there will be a certain concentration of salt.  To simulte this, we need to add ions.  We want to add enough ions to simulate a physiologically relevant concentration, and then add some extra ions to make sure our simulation has zero total charge.

We're going to use `genion` to add ions. Look up `genion` in the gromacs manual, as before, and see what it does, and what the flags mean.  

`genion` requires a `.tpr` file, so first we have to use `grompp` to make one.

~~~s
gmx grompp -f minimise.mdp -c 10_GlyT2_POPC_CLR_water_EM.gro -p GlyT2_POPC_CLR.top \
           -o 11_GlyT2_POPC_CLR_water_ions.tpr -maxwarn 1

gmx genion -s 11_GlyT2_POPC_CLR_water_ions.tpr -p GlyT2_POPC_CLR.top -conc 0.15    \
           -neutral -o 11_GlyT2_POPC_CLR_water_ions.gro
~~~

Choose a group that corresponds to the water molecules.

Gromacs will now replace water molecules with ions and update your `.top` file for you.

Do an energy minimzation of this new system.

~~~s
gmx grompp -f minimise.mdp -c 11_GlyT2_POPC_CLR_water_ions.gro -p GlyT2_POPC_CLR.top \
           -o 12_GlyT2_POPC_CLR_water_ions_EM1.tpr -maxwarn 1
    
gmx mdrun -v -deffnm 12_GlyT2_POPC_CLR_water_ions_EM1
~~~

Check the forces and energies.  Are they sensible?  If so, run a second EM.  If everything has gone well, this second EM should be quite short.

~~~s
gmx grompp -f minimise.mdp -c 12_GlyT2_POPC_CLR_water_ions_EM1.gro                 \
           -p GlyT2_POPC_CLR.top -o 13_GlyT2_POPC_CLR_water_ions_EM2.tpr -maxwarn 1

gmx mdrun -v -deffnm 13_GlyT2_POPC_CLR_water_ions_EM2
~~~

Check the forces and energies.  Have a look at it in vmd.  Does everything look correct?

If so, you are ready to do equilibration.

### Getting ready for equilibration

Equilibration allows our system to reach an equilibrium state, smoothing out artefacts from system preparation and reaching something that is more physically realistic.  If you look in vmd, there are some big artefacts in our system!  There's some big areas vacuum sitting between our protein and our membrane, and between our membrane and our solvent.  Furthermore, our solvent is highly ordered in a very unrealistic way.

To equilibriate our system, we are going to do a series of simulations where the protein is held in place, and everything else can relax around it.  We will slowly decrease the strength of the force constants holding the protein in place.

When we simulate a system, we control the temperature with a thermostat.  We usually control the temperature of the solvent separately to the temperature of the solute.  To do this, need to make an *index file* to tell groamcs which atoms are solvent and which are solute.

~~~s
gmx make_ndx -f 13_GlyT2_POPC_CLR_water_ions_EM2.gro -o GlyT2_POPC_CLR.ndx
~~~

This shows you all the different index groups that are already defined.  There should be one called `Water_and_ions`

We want to make a second group, which contains everything that is not in `Water_and_ions`.  To do this, we will use the `not` command.  

If `Water_and_ions` is group number 20, we will type

~~~s
    !20
~~~

and hit enter.  This will make a new group that contains every single atom that is not in the group `Water_and_ions`

hit enter again to see a list of all groups.There should now be a new group called `!Water_and_ions`

We need to change the name of that new group to `non_water_ions`.  If `!Water_and_ions` is group number 21, type

~~~s
    name 21 non_water_ions
~~~

and hit enter twice (once to save the new group, and a secon time to display the full list of groups)

If everything looks right, type `q` and then hit enter, and gromacs will save your new index file.

Next we need to edit `posre.itp`

Open `posre.itp` in your editor of choice.  It contains a series of restraints that apply forces to every protein atom to prevent it from moving away from a set location in the x, y and z dimension.  By default, that force is defined by the constant 1000.  Have a look in the manual to understand how position restraints work, and what kind of units they use.  We want to change the force constant for our position restraints from an explicit value (1000) to a variable we can change whenever we want, called `p_rest`.

Do a find and replace to change ` 1000  1000  1000 ` to ` p_rest p_rest p_rest `.  If you look in the different `equil_memb_N.mdp` files, you can see a statement that defines the value of `p_rest` for each simulation.  Think for a while about why we replaced ` 1000  1000  1000 ` with ` p_rest p_rest p_rest `, rather than just replacing ` 1000 ` with ` p_rest `.  There is a very important reason!

### Running equilibration

You are finally ready to run your equilibration!

Once again:  We want to series of 1 ns equilibration simulations, with position restraints starting at 1000, and decreasing to 10.  There should be five steps in total, in this order: 1000, 500, 100, 50, 10

Look at the file ` run_equil_1000_to_10.sh `

This is a script that will run the first two steps of your equilibration, the 1000 step and the 500 step.

How is the `grompp` command in this script different to the ones you have used previously?  What new flags are you using?  What do they mean?

You need to edit this file to add the remaining steps, for 100, 50 and 10.  You can use the 500 step as a template.  Think about why the 1000 step would not make a good template.

Once you have done that, you need to tell your computer that it is allowed to execute the program `run_equil_1000_to_10.sh` .  To do that, use the command `chmod`

~~~s
chmod +x run_equil_1000_to_10.sh
~~~

You might need to run that command as `sudo`

You are finally ready to run your equilibration.  This process will take several hours, so go through everything first and make sure it's all correct.  If you want, you can make a copy of your tutorial folder and run a test equilibration, where you edit the `equil_memb_N.mdp` files so they are very short by changing `nsteps` to 1000.  This will let you identify any errors in your `run_equil_1000_to_10.sh` script.

Once you are satisfied, run your eq with the command

~~~s
./run_equil_1000_to_10.sh
~~~

Good luck!  If all goes well, your final, fully equilibrated system, will be ready for you in a few hours.

While you are waiting, look at the gromacs manual entry for `pdb2gmx`.  What did it do?  What do the `-ter`, `-glu`, `-ss`, `-ignh`, `-vsite` `-heavyh` and `-deuterate` flags do?  Think about when you might want to use them.  Don't just assume the choices in this tutorial will be suitable for your own projects!

If everything worked, show your readme file and your finished system to the person administering the tutorial.  Explain to them what you did, and what problems you had on the way.  Pick one gromacs command you used, and explain how it works, and why it was useful.
