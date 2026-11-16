Contact information: taneli.tolppanen@gmail.com

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