# TwoSinglet_H3H2H1
Set up the asymmetric decays of a heavy scalar in Herwig 7

The ```unweighted_events.lhe.gz``` file was generated via MadGraph 5 using the TwoSinglet model (https://gitlab.com/apapaefs/twosinglet). 

```iota0``` is the heaviest scalar (925 GeV in the LHE file), which can then decay to a ```eta0``` of a choosen mass and a Higgs boson (h0), via:

$\iota_0 \rightarrow \eta_0 h$

The ```eta0``` particle was then set up to decay to W+W-.

See input file for details on the setup of the decay modes. 

Note that this is somewhat of a hack, and the built-in UFO interface is recommended! Works with Herwig 7.3.0.

## Herwig Model:

To run the MENewScalars model, place the ```NewScalars``` folder in the ```Herwig-7.3.0/Contrib``` directory. In that directory, type ```make``` to create the Makefile for the ```NewScalars``` model. 

Then, in the ```NewScalars``` directory, type ```make```.

An example input file exists: ```LHC-NewScalars.in```. 

IMPORTANT: Note that for now, the MHatMin and MHatMax cuts have to be set to be slightly below and slightly above the chosen mass for the ```iota0``` particle, otherwise this wil l note work!
