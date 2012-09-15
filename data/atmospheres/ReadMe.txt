J/A+A/333/231       O-M stars model atmospheres (Bessell+ 1998)
================================================================================
Model atmospheres broad-band colors, bolometric corrections and temperature
calibrations for O - M stars
       Bessell M.S., Castelli F., Plez B.
      <Astron. Astrophys. 333, 231 (1998)>
      =1998A&A...333..231B      (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Models, atmosphere ; Photometry, UBVRIJKLMNH
Keywords: stars: atmospheres; fundamental parameters; general

Abstract:
    Broad band colors and bolometric corrections in the
    Johnson-Cousins-Glass system (Bessell, 1990PASP..102.1181B; Bessell &
    Brett, 1988PASP..100.1134B) have been computed from synthetic spectra
    from new model atmospheres of Kurucz (1995a, priv. comm.), Castelli
    (1997, priv. comm.), Plez, Brett & Nordlund (1992A&A...256..551P),
    Plez (1995-97, priv. comm.), and Brett (1995A&A...295..736B,
    1995A&AS..109..263B). These atmospheres are representative of larger
    grids that are currently being completed. We discuss differences
    between the different grids and compare theoretical color-temperature
    relations and the fundamental color temperature relations derived
    from: (a) the infrared-flux method (IRFM) for A-K stars (Blackwell &
    Lynas-Gray 1994A&A...282..899B; Alonso et al. 1996A&AS..117..227A) and
    M dwarfs (Tsuji et al. 1996A&A...305L...1T); (b) lunar occultations
    (Ridgway et al. 1980ApJ...235..126R) and (c) Michelson interferometry
    (Di Benedetto & Rabbia 1987A&A...188..114D; Dyck et al.
    1996AJ....111.1705D; Perrin et al. 1998A&A...331..619P) for K-M
    giants, and (d) eclipsing binaries for M dwarfs. We also compare
    color-color relations and color-bolometric correction relations and
    find good agreement except for a few colors. The more realistic fluxes
    and spectra of the new model grids should enable accurate population
    synthesis models to be derived and permit the ready calibration of
    non-standard photometric passbands. As well, the theoretical
    bolometric corrections and temperature-color relations will permit
    reliable transformation from observed color magnitude diagrams to
    theoretical HR diagrams.

File Summary:
--------------------------------------------------------------------------------
  FileName   Lrecl Records Explanations
--------------------------------------------------------------------------------
ReadMe          80       . This file
table1.dat      82     224 Color indices and bolometric corrections from ATLAS9
                             overshoot models (Kurucz, 1995a, priv. comm.)
table2.dat      81     224 Color indices and bolometric corrections from ATLAS9
                             no-overshoot models (Castelli, 1997, priv. comm.)
table3.dat      82     252 Color indices and bolometric corrections from ATLAS9
                             models for Teff>=9000K (Kurucz, 1993)
table4.dat      71      93 Color indices and bolometric corrections from NMARCS
                             giant branch models of Plez, Brett and Nordlund
                              (1992A&A...256..551P)
table5.dat      74     186 Color indices and bolometric corrections from NMARCS
                             giant branch models of Plez (1995, priv. comm.).
table6.dat      74      76 Color indices and bolometric corrections from NMARCS
                             M dwarf models of Plez (1997, priv. comm.)* and
                             Brett (1995A&A...295..736B and 1995A&AS..109..263B)
--------------------------------------------------------------------------------

See also:
     VI/39 :   Model Atmospheres (Kurucz, 1979)

Byte-by-byte Description of file: table1.dat table2.dat table3.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  5  I5    K       Teff      Effective temperature
   7- 10  F4.2 [cm/s2]  logg      Surface gravity
  12- 16  F5.2  mag     BC(K)     Bolometric correction in Kmag
  18- 22  F5.2  mag     BC(V)     Bolometric correction in Vmag
  24- 29  F6.3  mag     U-B       U-B colour index
  31- 36  F6.3  mag     B-V       B-V colour index
  38- 43  F6.3  mag     V-R       V-R colour index
  45- 50  F6.3  mag     V-I       V-I colour index
  52- 57  F6.3  mag     V-K       V-K colour index
  59- 64  F6.3  mag     J-H       J-H colour index
  66- 71  F6.3  mag     J-K       J-K colour index
  73- 78  F6.3  mag     K-L       K-L colour index
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table4.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  4  I4    K       Teff      Effective temperature
   6- 10  F5.2 [cm/s2]  logg      Surface gravity
  12- 15  F4.2  solMass Mass      Mass
  17- 20  F4.2  mag     BC(K)     Bolometric correction in Kmag
  22- 26  F5.2  mag     BC(V)     Bolometric correction in Vmag
  28- 32  F5.3  mag     B-V       B-V colour index
  34- 38  F5.3  mag     V-R       V-R colour index
  40- 44  F5.3  mag     V-I       V-I colour index
  46- 50  F5.3  mag     V-K       V-K colour index
  52- 56  F5.3  mag     J-H       J-H colour index
  58- 62  F5.3  mag     J-K       J-K colour index
  64- 68  F5.3  mag     K-L       K-L colour index
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table5.dat table6.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label  Explanations
--------------------------------------------------------------------------------
   1-  4  I4    K       Teff   Effective temperature
       5  A1    ---     Note   [*] *: Plez 1997 (priv. comm.) data in table6.dat
   6- 10  F5.2 [cm/s2]  logg   Surface gravity
  12- 16  F5.2  Sun     [Fe/H] Metallicity
  18- 21  F4.2  mag     BC(K)  Bolometric correction in Kmag
  23- 27  F5.2  mag     BC(V)  Bolometric correction in Vmag
  29- 33  F5.3  mag     B-V    ? B-V colour index
  35- 39  F5.3  mag     V-R    V-R colour index
  41- 45  F5.3  mag     V-I    V-I colour index
  47- 52  F6.3  mag     V-K    V-K colour index
  54- 58  F5.3  mag     J-H    J-H colour index
  60- 64  F5.3  mag     J-K    J-K colour index
  66- 70  F5.3  mag     K-L    K-L colour index
--------------------------------------------------------------------------------

Acknowledgements: Fiorella Castelli <castelli@astrts.oat.ts.astro.it>

References:
   Kurucz, R.L. 1993, "ATLAS9 Stellar Atmosphere Programs and 2 km/s grid",
                        CD-ROM No 13
================================================================================
(End)                                         Patricia Bauer [CDS]   05-Mar-1998
