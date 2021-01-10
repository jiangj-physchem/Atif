## Atif
Atif is an advanced theoretical tool for inhomogeneous fluids. Atif can be used to study many kinds of fluids near an interface:
* hard sphere fluids with size symmetry/asymmetry 
* electrolytes including size symmetry/asymmetry and monovalent/multivalent ions
* flexible and semiflexible uncharged/charged polymers (polyelectrolytes)
* block uncharged/charged polymers and sequential uncharged/charged polymers

Currently, this theoretical tool is only for open systems (Grand-canonical ensemble). Howerver, it is very easy to make use of this tool to study issues in a closed system (canonical ensemble) with little changes (please directly contact [the developer](https://github.com/jiangj-physchem) if necessary). In addition, the current version can only study the inhomogeneous properties near a **planar surface**. 

In the future, we would also like to involve **spherical** and **cylindrical interfaces**.

## The theories we used in this tool

**I. self-consistent field theory (SCFT)** 
* the flexible polymer is modeled by freely jointed chain model
* the semiflexible polymer is modeled using discrete worm-like chain model
* the excluded volume effect is treated using local incompressible condition
* the electrosatic potenital is obtained using point charge model.

**II. density functional theory (DFT)**
* the flexible polymer is modeled by freely jointed chain model
* the semiflexible polymer is modeled using discrete worm-like chain model;
* the excluded volume effect is involved using modified fundamental measure theory (MFMT)
* the electrostatic potential is obtained using truncated shell model (TSM)
* the electrostatic correlations are involved using a functional mean spherical approximation (MSA)
* the non-bonded chain connectivtiy contirbutions are considered using thermodynamic perturbation theory (TPT)

I will cite the corresponding references subsequently.

## How to use Atif?
Before you enjoy in **Atif**, please read the [GETTINGSART.md](GETTINGSART.md) file 

## Want to thank the developer?
If you think Atif really help you during your research career, please feel free to thank the developer in the "Acknowledge" section in your publications. If you would like to put the developer in the author list of your publications, please directly contact [the developer](https://github.com/jiangj-physchem).

## Contributing

Currently, we have not set up the rules for submitting pull requests. Please directly contact [the developer](https://github.com/jiangj-physchem) if you would like to contribute.

## License

Atif is licensed under the [MIT License](LICENSE.md).
