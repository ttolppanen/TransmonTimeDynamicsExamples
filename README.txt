Contact information: taneli.tolppanen@gmail.com
Info about the code I used in my projects is at the bottom.

These codes use my custom packages, which can be installed from my registry.
To add the custom registry, in julia REPL use ] to go to the package manager mode (PKG), and type 

registry add https://github.com/ttolppanen/ttolppanen_JuliaRegistry.git

To activate this projects environment, all in PKG mode type 

activate . (it finds the environment from the folder the REPL is active in)
update
instantiate

Then you should have all packages installed and ready to run the codes.

Some relevant packages available from my registry are:
QuantumStates (https://github.com/ttolppanen/QuantumStates)
QuantumOperators (https://github.com/ttolppanen/QuantumOperators)
QuantumTimeEvolution (https://github.com/ttolppanen/QuantumTimeEvolution)

There can be more usefull things inside these packages that come up 
here in the examples. You can look whats in the packages by following 
the link, and going to the "src" foulder. In the "test" foulder, there 
can be more usefull examples and how to use the packages.

There is also:
PlotAndSave (https://github.com/ttolppanen/PlotAndSave)

which saves the data along with making a simple plot for it, I found it
very usefull since you often need to save the data while making the plots.

NOTE THAT I USED JULIA VERSION 1.11.6 (or older) WHEN MAKING THESE PACKAGES. They should 
work fine, but if you run into issues, consider installing an older version of julia.
Also, I strongly recommend using Visual Studio Code.

I want to to highlight that I feel like the implementation here for calculating the Bose-Hubbard Hamiltonian
is efficient, so consider using that or taking some inspiration if you run into problems with scaling that for
large systems in your implementation. The main trick that I found useful was to use the sparse matrix contructors.
This together with Olli's implemenation where you only loop through the basis states once.

---ABOUT PROJECT CODE---
More code related to Paper I measurement induced phase transitions
https://github.com/ttolppanen/Measurement-Induced-Phase-Transition-in-Transmons

For Paper II transmon array based leakage reduction unit (LRU)
https://github.com/ttolppanen/TestingStuff
Here the relevant code can be found from
TestingStuff/Physics/LeakageRemovalUnit/PlotsForPaper

For Paper III leakage in tunable couplers
https://github.com/ttolppanen/LRUFollowUpPaper
Here the relevant code can be found from
LRUFollowUpPaper/LocalPackages
LRUFollowUpPaper/LRUFUP_Vault/Data/0PAPER
and in 0PAPER:
LeakageHopping200 and LeakageHopping201 relate to paper plots Fig.3 a-d and e-h respectively.
LeakageTunneling to Fig. 4.
LRUComparison has code for simulating the LRUs we studied in the paper. 
