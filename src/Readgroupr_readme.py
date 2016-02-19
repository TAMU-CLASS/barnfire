'''
Andrew Till
Summer 2014

Read in a GROUPR file from NJOY and store in a useful data structure.

--------------------------------------
GROUPR stores reactions by (MF, MT) pairs.
MF is the type of value to read
MF == 3 --> cross section (vector)
MF == 6 --> transfer matrix (matrix)

MT is the type of reaction to read
MT == 1 --> total
MT == 2 --> elastic scattering
MT == 18 --> fission (no nu)
MT == 452 --> total nu
MT == 221 --> free gas thermal
etc.

--------------------------------------
The GENDF output used by GROUPR is both good and bad. On the plus side, it is very modular, with each reaction being essentially independent. On the minus side, it does a poor job of giving sizes of arrays upfront, especially for matrices (which are stored somewhat sparsely), and number of temperatures. In the latter case, you basically have to read the whole file to know which reactions are present at which temperatures before you know how much memory to allocate.

GROUPR has many quirks. One of these quirks involves asking for the transfer matrix for the (n,2n), (n,3n), (n,4n), and continuous inelastic scattering. Although the reaction is in MF 6 and the resulting location in the GENDF file is in MF 6, the proper input to GROUPR is to ask for the reactions with MF 8. This is not an issue if all reactions are desired ('3/' and '6/').

Another quirk of GROUPR is that it stores its energies from low to high, whereas most multigroup cross sections are stored high to low. The code switches order internally after it reads the cross sections and before saving them.

Yet another quirk of GROUPR is that the first background cross section ($\sigma_0$) value must be infinity. Subsequent values must be in descending order.

--------------------------------------
NJOY uses the concept of background cross section to treat within-group spectral shielding effects in the resolved and unresolved resonance regions. The spectrum for nuclide $i$ is assumed to be $C(E) / (\sigma_{t,i}(E) + \sigma_{0,i})$, where $C(E) = 1/E$ usually. This is valid under the narrow resonance (NR) approximation. The value of background cross section, $\sigma_{0,i}$ comes from the effect of the other nuclides in the material and, optionally, an escape cross section from equivalence theory. In practice, NJOY is run with multiple values of the background cross section. Reactions are then interpolated on the $\sigma_0$ grid for each material composition for each group ($\sigma_{0,i} becomes \sigma{0,i,g}$). This process requires iteration, often referred to as Bondarenko iteration. More advanced methods will reevaluate the escape cross section (say, via independent fixed-source problems for each group) at each stage of this iteration, though that is beyond the scope of this work.

The background cross section only affects to resolved and unresolved resonance regions. Background cross sections are meant to capture the effect of resonances not resolved by the group structure. At thermal and epithermal energies, there are either no energies, or they are resolved by the group structure (wide resonances). At fast energies and energies near the upper portion of the unresolved resonance region, there are either no resonances, or the resonances are so small and finely spaced that they act as a continuum.

For the above reasons, not all reactions are stored at multiple background cross section values, only those which have non-zero cross section in the resolved (or lower portion of the unresolved) resonance range. This means threshold reactions, such as inelastic scattering and (n,2n), are only stored for the infinite-dilute case. Similarly, thermal scattering is stored at only one backround cross section. While the fission cross section may be self-shielded, the spectrum is not. For low energies, the reaction is dyadic (constant spectrum) which is rigorously independent of how the within-group spectrum is handled. For high energies, although the spectrum is dependent on the incident group, these energies are above the unresolved resonance range and the fission cross section here is independent of background cross section.

The reactions stored at multiple background cross sections are
Name              (MF, MT)
total             (3, 1)
elastic           (3, 2) and (3, 6)
fission           (3, 18)
radiative capture (3, 102)

Treatment of the unresolved resonance range is handled through the use of the background cross section. In this range, the actual positions of the cross sections themselves is not known, only the distribution of cross sections' positions, widths, and magnitudes. In this range, the narrow resonance approximation is valid, and the use of background cross sections will produce an accurate cross section. This is good news, as the PG-FEMG method requires partitioning the energy domain, which requies knowing where the cross sections are.

--------------------------------------
Not all reactions need to be stored at multiple temperatures. Temperatures only affect resonances and thermal treatment. As resonances appear only in the resolved and unresolved resonance regions, the same arguments used above apply to which reactions are affected by changing temperatures, modulo the thermal cross sections.

GROUPR treats multiple temperatures more loosely than multiple background cross sections. The user is free to specify whichever reactions at whichever temperatures they wish. Including reactions which are unaffected by increases in temperature is a poor use of resources, both with respect to computing and storage.

The reactions which ought to be stored at multiple background cross sections are
Name              (MF, MT)
total             (3, 1)
elastic           (3, 2) and (6, 2)
fission           (3, 18)
radiative capture (3, 102)
Thermal XS        (3, 221)-(3,246) and (6, 221)-(6,246)

--------------------------------------
Basic resonance interference effects are treated by using the flux calculator capability in NJOY coupled to the admixed moderator method with U-238. This latter method assumes U-238 is the dominant resonant nuclide, which is valid for most thermal reactors.

--------------------------------------
--------------------------------------
GROUPR treats cross sections with non-unity multiplicity in a consistent, if not intuitive, way. The cross section (MF 3) is calculated without the multiplicity. The transfer matrix (MF 6) is calculated with the multiplicity. The total cross section remains the sum of the partial cross sections, while the transfer matrix can be used as-is, without the need to tack on factors afterward. For example, (MF, MT) = (3, 16) ((n,2n) cross section) has no factor of 2, but (6, 16) ((n,2n) transfer matrix) has a factor of 2. Sums over outgoing groups (columns) of (6,16) equals 2 times (3, 16) for the same incident group. Similarly, (3, 18) (fission cross section) has no factor of nu, though it can be asked for explictly using (3, 452). (6, 18) (fission transfer matrix) has a factor of nu. This presents difficulties if for uninformed converter codes.

All transfer matrices are given in Legendre moments of the cosine of the scattering angle. Weighting the higher-order moments is NOT done with the same flux as is used for everything else, but rather $C(E) / (\sigma_{t,i}(E) + \sigma_{0})^(l+1)$, where $l$ is the Legendre moment. All transfer matrices except fission, which is isotropic, are given to the user-desired number of moments.

The thermal scattering cross sections do not necessarily have the same magnitude as the elastic scattering cross sections. Since the thermal XS are meant to overwrite the elastic XS, the total XS should be updated to be consistent with this new source. For many materials, there is an elastic and inelastic component to the thermal XS. These should be added together before they are used to overwrite the elasitc cross sections. For the transfer matrices, at low energies, the original contributions include only the elastic cross section, so the thermal cross sections can be substituted without worrying about missing something.

--------------------------------------
The thermal scattering treatment in NJOY is done with the THERMR module. Each thermal type is processed independently and added to the PENDF tape. A maximum energy for the thermal treatment is specified for each THERMR run so long as the desired temperature is below 3000 K. The sparsity structure of the thermal scattering transfer matrix is in general _different_ for each temperature, which requires a pre-loop to go through and figure out sizes for allocation purposes.

Also, need to handle some reactions not showing up for some T (thermals).

--------------------------------------
This module internally stores matrices in a CSC-like (compressed storage column) format. Transfer matrices are banded by column and generally full within the skyline. The values for column (which corresponds to one incident group) are stored contiguously. A pointer index array is used to determine where in the values array the cross sections are located. An index array provides the lowest-index group to which scattering occurs. This last point differs from the standard CSC format but takes advantage of the fullness within the skyline. During the read, we must also store the incident group, because certain reactions, such as threshold reactions or thermal reactions, are not populated for all groups.

The fission transfer matrix (6, 18) is stored in a compressed format in GENDF. This format consists of three sections. In the first section, the fission spectrum is stored for the low-energy groups. In the second section, the production cross section (nu $\sigma_f$) is stored for the low-energy groups. In the third, the fission matrix is stored for the high-energy groups.

Internally, I store the fission matrix in one piece: the low-energy spectrum in the first column, and the fission source (now normalized to nu $\sigma_f$ instead of 1) in the remaining columns. As the fission spectrum has no $\sigma_0$ or temperature dependence, I store it in a 2D array (g, g') (numpy uses C-like ordering). Knowing the size of the matrix (g) and the total number of groups, I can trivially calculate which groups use their own spectrum.

Most transport codes, including PDT, cannot read hybrid fission matrices. They instead use a fixed fission spectrum for the entire incident energy range. I compute an effective chi as

\chi_{g,eff} = ( \chi_{g,low} \sum_{g' \in low} nu \sigma_{f,g'} \phi_{g'} +
                              \sum_{g' \in high} nu \sigma_{f,g' -> g} \phi_{g'} )  /
                \sum_{g'} nu \sigma_{f,g'} \phi_{g'}

This requires reading in (3, 18), as I choose not to store nu $\sigma_f$ from (6, 18).

--------------------------------------
The reason NJOY doesn't work well for fast reactors is that it doesn't do self-shielding all that well, especially for resonance interference effects. You can get basic resonance interference effects by considering a U-238 background. However, NJOY is more atomistic than the material level, so it doesn't know about mixing, which means you can't get those effects correct. MC^2-3 does a much better job, but may be more expensive.

get_reaction_list (with modification) works for ENDF files ;)

--------------------------------------
NJOY does a somewhat poor job of handling many-body (three or more) threshold reactions (negative Q value) when it comes to producing a multigroup transfer matrix. It all stems from the ENDF data representation, which often gives the secondary energy distributions for a discrete number of incident energies, one of which is the effective threshold energy. For incident particles at this energy, the ejectiles have almost no energy, as all of it is soaked up in forming the product nucleus. The secondary energy distribution for such incident energies particles is very low. Often, the next-highest incident energy data point is sufficiently higher that its secondary energy distribution is orders of magnitude higher than for the previous case. Further, the lower cutoff for this case is often _above_ the upper cutoff of the previous case. When a fine group structure is used, this leads to a _gap_ in the groupwise secondary energy distribution in the transfer matrix, which is a nonphysical artifact of NJOY using a low-fidelity data representation for its threshold reactions. What is physical is that threshold reactions can contribute non-negligibly to low secondary energies for incident energies near the effective threshold values.

Notice this is not a problem for reactions that do not have thresholds (negative Q values) and that this is not a problem for 2-body reactions. 2-body reactions have a nice analytic representation of the exiting energy as a function of the angle of scattering and are stored differently. Finally, this issue of storing the secondary energy distribution on a lower energy grid for the lowest-energy incident particle seems only to apply to MT 16 and 17 (n,2n) and (n,3n) for minor actinides. MT 91, (n,n_c), does not have this problem because the lower bound of all secondary energy distributions is the same for all incident energies in the evaluations investigated.


'''



#TODO: Finish creating xs readers, esp. fission

