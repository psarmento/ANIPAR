
** DISCARDED because the original source for those abundances was not found or the values do not match:

  SPECTRUM/stdatom.dat:
	Solar atomic abundances from Grevesse & Sauval (1998) [1] provided by SPECTRUM[2]. 

  Kurucz/stdatom.dat:
	Solar atomic abundances from Kurucz atmosphere models (February 2012).

** REST:

Castelli/stdatom.dat:
    Solar atomic abundances from Castelli atmosphere models (August 2012).

MARCS/stdatom.dat:
    Solar atomic abundances from MARCS atmosphere models (September 2012).
    Transformation: value = value[MARCS] - 12 - 0.04


""" Contents (from SPECTRUM manual):
The first column is the atomic or molecular code, which gives both the atomic number and the ionization state. The code ``26.0'' refers to an Fe I (neutral iron) line, as the atomic number of iron is 26, and the neutral ionization state is given as decimal 0. The code for Fe II is 26.1, the code for Ca III is 20.2 and so on. For diatomic molecules, the code specifies the atomic makeup of the molecule. Thus, H is 101.0, the two ``1''s referring to the two hydrogens, CH is 106.0, CO 608.0, MgH 112.0, TiO 822.0, etc. The lightest element always comes first in the code, so that 608.0 cannot be confused with NdO, which would be written 860.0.

The second column contains, for the atoms, the logarithmic solar abundance in terms of number densities relative to the total number density: log(A/Ntotal). This is the way that the abundance scale is defined in ATLAS9 and 12, but it is not the standard way of representing elemental abundances (expressed with respect to hydrogen, on a scale in which the log of the abundance of hydrogen is set equal to 12.0).

	** SPECTRUM uses a scale in which the abundances are expressed with respect to the total number density of atoms (and ions). To convert between the two scales is easy; use the following formula

		log(A/NH) = 12.0 + log(A/Ntotal) - log(NH/Ntotal)
        log(A/Ntotal) = log(A/NH) - 12.0 + log(NH/Ntotal)

	For log(NH/Ntotal), SPECTRUM uses the hydrogen abundance in the atomic and molecular data file (stdatom.dat or a corresponding custom file). The value in the current version of stdatom.dat is -0.0360.

	** For the molecules, this column is set equal to -7.50 for obscure reasons.

The third column contains the mean atomic mass for both atoms and molecules. 

The fourth column contains for atoms the first ionization energy in electron volts (eV); for molecules this column contains the dissociation energy, also in eV.

The fifth column contains for atoms the second ionization energy and for molecules the "reduced mass" ((m1 * m2) / (m1 + m2 ). 

The sixth column contains, for the atoms, the third ionization energy, and for the molecules a "fudge" factor by which the "gf" value for each spectral line is multiplied. This gives the user a very rough and blunt instrument for modifying molecular band strengths, for instance, after their appearance in the solar spectrum.

The seventh column gives, for the atoms the fourth ionization energy and is identically zero for the molecules. 

Finally, the eighth column gives the maximum charge (the highest ionization state) supported by SPECTRUM. For instance, a "1" means that SPECTRUM can compute spectral lines only for neutral and singly ionized species (I and II), a "3" means that SPECTRUM has support for ionization states I, II, III and IV. For molecules this entry is 0. All naturally occurring elements are supported (hydrogen - uranium) by SPECTRUM and 15 diatomic molecules. Adding support for another molecule is not as simple as adding another line to stdatom.dat. It also involves code changes, addition of the molecular partition function, etc.
"""

SPECTRUM will use the value for [M/H] to scale the abundances of the metals (everything except for H and He) in the atomic and molecular data file stdatom.dat.

Thus, if [M/H] = -1.0, SPECTRUM subtracts 1.0 from the logarithmic abundances in stdatom.dat for lithium through uranium. Since the opacity due to the iron peak elements is so important in establishing the structure of the stellar atmosphere, it is advisable to use a stellar atmosphere model with [M/H] as close as possible to the [Fe/H] of the star you are analyzing.

For instance, in the sun (and in stdatom.dat), the logarithmic abundance of calcium (code = 20) is -5.68. Let us suppose that you actually want to use an abundance of -6.38. If the [M/H] of the model is , then the [M/H] scaling changes the abundance from stdatom.dat to -6.68, so this means that to get an abundance of -6.38 you will need to adjust the abundance in stdatom.dat to -5.38.


[1] Grevesse, N. & Sauval, A.J. 1998 Space Science Reviews 85, 161 
[2] http://www1.appstate.edu/dept/physics/spectrum/spectrum.html
