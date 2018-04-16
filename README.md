# iondosecvtr

A simple Python utility to convert relative depth ionisation to relative depth dose, in line with the formalism outlined in the Code of Practice. The empirical formula used to compute the water-air mass stopping power ratios is valid for clinical electron beams with QI between 1 cm and 20 cm and for water depths between 0.02×R<sub>50,D</sub> and 1.2×R<sub>50,D</sub>.

# usage

The conversion of depth ionisation to depth dose may be performed by *e.g.*

```Python
>>> iondosecvtr.cvtiondose(file='6_MeV.csv',r50d=2.34)
```
where the input CSV file must be formatted as shown in [6_MeV.csv](https://github.com/archon88/iondosecvtr/blob/master/data_files/6_MeV.csv); notably, the depth is assumed to be expressed in millimetres, while QI (r50d) is in centimetres, following the convention used in the Code of Practice. This function outputs a second CSV file with two additional columns containing the water-air mass stopping power ratio (SPR) for the water depth and electron beam quality, and the relative dose at depth. Both the relative ionisation and relative dose are expressed as percentages, normalised to their respective maxima.

Then, plots of the relative dose, relative ionisation, and SPR against depth may be created by passing the name of the file output by **cvtiondose** and the nominal beam energy (in MeV):

```Python
>>> iondosecvtr.makeplots(file='6_MeV_withdose.csv',energy=6)
```
