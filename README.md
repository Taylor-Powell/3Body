# Three-Body Integral Equations (Refactored)

## Purpose
This code is intended to be a compilation of the various programs I wrote during research with Dr. Raul Briceno and Dr. Andrew Jackura during undergrad. I will list functionality sorted by executables located in _main/_

### free_2part.cpp
This is the simplest code in use. A sample input parameter file is located in _inputs/params_free2part.dat_.

The code instantiates the Free2Part::Data class, which takes the parameter file as an input. The class is defined in _Free2Part.h_. The constructor sets a few private variables and then reads the parameters from the input file, throwing an error if any necessary values are omitted. The constructor also sets the output filename based on provided parameters (this avoids potential name conflicts from generic naming of output files).

After creating the class, the main function calls the _free_2_out()_ function. This iterates over each box side length $L$ and calls _gen_free_spec_2part()_ from _FreeSpec.h_. The non-interacting spectrum for each box size for the particle masses provided is then computed according to

$$E_n=\sqrt{m_1^2+\left(\frac{2\pi}{L}\vec{n}_1\right)^2}+\sqrt{m_2^2+\left(\frac{2\pi}{L}\vec{n}_2\right)^2}$$

The total energy of the system is bounded by $E_{\text{max}}$. Since we are only interested in the low-energy spectrum, a (generous) cutoff in $\vec{n}=(n_x,n_y,n_z)\in\mathbb{Z}^3$ is included, where each component is restricted as $|n_i|\le5$. 

After the spectrum is computed for each volume, it is sorted and printed to the output file for later printing. As an example, formatted lines of the output file for the default parameters provided are:

| $L$ | $E_1$ | $E_2$ | $E_3$ | $E_4$ | $E_5$ | $E_6$ | $E_7$ | $E_8$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 2 | 2.54311 | 2.98911 | 3.37671 | 3.72419 | 4.04191 | 4.3364 | 4.87229 |
| 7 | 2 | 2.68751 | 3.23194 | 3.69705 | 4.10986 | 4.48482 | 4.83077 | | 
| 6 | 2 | 2.89594 | 3.57393 | 4.1424 | 4.64176 | | | 

As seen in the table, as the box size decreases, the number of accessible energy levels below $E_{\text{max}}$ will necessarily decrease. The output is formatted so the first row contains the most energy levels, which makes it easier to plot in many scripting languages.

### free_coupled22.cpp
This runs in a similar manner to _free_2part.cpp_, except now the input parameter file is split into global, volume, and particle parameters. 

After initializing the Coupled22::Data class, the constructor again reads the input file and throws an error if something wasn't properly initialized. After creating the class, the main function calls the _freeSpecOut()_ function, which builds and outputs the non-interacting spectrum for each box volume to two separate output .dat files in _data/coupled22/_ for later use.

### intspec_coupled22.cpp
Work in progress.
