
       solvRadialSchroedinger
          
     C++ programme
     Very simple code using Object-Oriented programming
------------------------------------------------------------
------------------------------------------------------------

 - Solution of the scattering problem in position space.
 - For spherical potential (Partial Waves Decomposition) 
 - Currently using only a Woods-Saxon potential.
 - In position space, the Schroedinger Equation is solved
   using the numerov algorithm.

-----------------------------------------------------------
 inputs:
      Ecm = Energy [MeV]
      l   = Orbital angular quantum number
   Calculates the scattering wave and uses it.
   To obtain the spectral function.

   Default stuff:
  Reduced mass is assumed to be the Neutron Mass.
  The Woods-Saxon potential is good for a nuclues 
  like 40Ca, therefore in practice the code only 
  solves  n+40Ca scattering.
    ... huge room for improvement
