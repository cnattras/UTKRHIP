// -*- C++ -*-
#include "Rivet/HeavyIonAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet
{
/// **********************************************************************************
/// This is a fake analysis calculating correlation functions off of trigger particles.
/// This code is written as if the reference paper has the following methodology:
///   It plots over events in the 10-30% and 50-80% centrality range, but that data is seperate
///   Triggers particles are defined as Pi0s with 13 GeV < pT < 20 GeV or
///      gamma particles with 8 GeV < pT < 20 GeV. Both are |eta| < 1.0
///   Associated partices are all charged hadrons with pT > 1.2 GeV and less than trigger pT


	class STAR_2014_I1234567 : public HeavyIonAnalysis
	{
	// Variables pertaining to just this class
	// I moved the private section above the public to introduce
	//   the variables you will see in init and later
	private:
		const int N_CENT_TYPES = 2; // number of centrality bins that you need to worry about, this case two, 10-30, 50-80
		const int CENT_TYPE_EDGES[N_CENT_TYPES][2] = {{10,30},{50,80}}; // 10-30, 50-80
		
		// histograms; Yours will almost certainly be different
		Histo1DPtr _h1dPhi[N_CENT_TYPES]; // delta phi, split by centrality
		
		// Number of trigger particles for each centrality over
		//  all events. Needed in finalize() for normalization
		unsigned long long int nTrigger[N_CENT_TYPES];

	public:
		/// Constructor
		STAR_2014_I1234567():HeavyIonAnalysis("STAR_2014_I1234567") {}

		/// Book histograms and initialise projections before the run
		///   Every histogram you will use must be booked, and every projection
		///   you will use must be declared
		void init()
		{
			//**** Select centrality method ****
			// The 50 here is the number of events used JUST to determine future events' 
			//   centrality. These events will show a centrality of -1.0, which is not
			//   physical. You can adjust this number for testing, but never have it lower
			//   than the number of events you are testing against.
			//   The test file ampt.hepmc has 100 events
			addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "IPMethod");

			//**** Trigger particle set ****
			// Many of you will have "trigger" particles. This is code to trigger off of a 
			//   specific type or particle in a certain range. There will ALWAYS be an "abseta"
			//   cut on all of your particles. It may differ from your trigger and associated 
			//   particles. This example uses Pi0 and gamma of differet energies as triggers.
			const int pidPI0 = 113; // The PID for Pi0
			const int pidGAMMA = 22; // The PID for Gamma
			// Find other pid codes at: http://home.fnal.gov/~mrenna/lutp0613man2/node44.html
			// I apply the following cuts for trigger particles:
			//  |eta| < 1.0  for all particles
			//  and it is either
			//      A Pi0 with pT between 13 and 20 GeV
			//  or  A gamma with pT between 8 and 20 GeV
			Cuts cutTrigger = Cuts::abseta < 1.0 && ((Cuts::pid == pidPI0 && Cuts::pt > 13.0 * GeV && Cuts::pt < 20.0 * GeV) || (Cuts::pid == pidGAMMA && Cuts::pt > 8.0 * GeV && Cuts::pt < 20.0 * GeV))
			FinalState fs(cutTrigger);
			declare(fs, "partTrigger");

			//**** Associated particle set ****
			// Charged Particles with cuts
			// Max pT for an associated particle will be limited by the pT of the trigger
			//  but that is based on the specific trigger. So that check will be done later
			//  inside a loop in analysis
			Cuts cutAssoc = Cuts::abseta < 1.0 && Cuts::pt > 1.2 * GeV && Cuts::pt < 20.0 * GeV;
			ChargedFinalState cfs(cutAssoc);
			declare(cfs, "partAssoc");				

			//**** Book Histograms ****
			//  This case I only have two since this not a real analysis
			//  I am looping over all centrality ranges I care about
			for (int i = 0; i < N_CENT_TYPES; ++i) {
				// This creates a histogram with the name "d01-x01-yNN" where NN = i
				//  and uses the reference data from the .yoda file associated with
				//  this analysis
				_h1dPhi[i] = bookHisto1D(1,1,i); 
			}
			
			//**** Initialize counters ****
			// Set number of events in counters to zero
			// These are used later for normalizing histograms
			//  in finalize
			for (int i = 0; i < N_CENT_TYPES; ++i) nTrigger[i] = 0;
		}

		/// Per event calculations
		///  Do your per event calculations and
		///  fill your histograms here
		void analyze(const Event& event)
		{
			// Get the centrality for each event
			const double c = centrality(event, "IPMethod");

			// The first 50 events (number from addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "IPMethod") )
			//  will give a centrality outside of 0-100, specifically -1.0
			if ((c < 0.) || (c > 100.)) vetoEvent; 
			
			// Alternatively, you may want to set the acceptance to a narrower range based on what your paper plots
			// For my fake paper there are two ranges. You can use these as indices for an array of histograms instead.
			int centralityIndex = -1; // -1 - not in plot range
									 // 0 - between 10-30
									 // 1 - between 50-80
			// Find which centrality range of interest this falls into
			for (int i = 0; i < N_CENT_TYPES; ++i) {
				if (c > CENT_TYPE_EDGES[i][0] && c <= CENT_TYPE_EDGES[i][1]) {
					centralityIndex = i;
					break; // since theres no centrality overlap, break the loop when you find the correct index
				}
			}
			// If not a centrality of interest, stop processing the event
			if (centralityIndex == -1) vetoEvent;
			
			

			// In the background RIVET is handling the the declared projections
			//  but you need to access them inside analyze to use tCFShem
			const FinalState & fs = apply<FinalState>(event,"partTrigger");
			const ChargedFinalState & cfs = apply<ChargedFinalState>(event,"partAssoc");
			
			// Now that you have the projection objects, you need to get the particles from them.
			Particles tracksTrigger = fs.particlesByPt();
			Particles tracksAssoc = cfs.particlesByPt();
			
			// Increase appropriate counter
			nTrigger[centralityIndex] += tracksTrigger.size();
			
			// foreach goes through each Particle in tracksTrigger, and references it in partTrigger.
			//  Its a simplification of a for loop where you would have had to find the number of entries
			//  and loop over them.
			// I am assuming I am calculating delta-phi between trigger particles and charged hadrons
			double deltaPhi;
			foreach (const Particle& partTrigger, tracksTrigger) {
				// Loop over all associated particles
				foreach (const Particle& partAssoc, tracksAssoc) {
					// Only include associated particles with pT less than the trigger
					if (partAssoc.pt() < partTrigger.pt()) {
						deltaPhi = partAssoc.phi() - partTrigger.phi();
						// In this fake paper plots run from 0 - 2 pi
						//   so make sure deltaPhi falls in this range
						// M_PI is part of the <math.h> header. You need to use the define: #define _USE_MATH_DEFINES
						//   See header includes at top of file for example
						while (deltaPhi < 0) deltaPhi += 2 * M_PI; 
						
						// Fill histogram with delta phi information
						//  In this case use a 1 as the second parameter, the weight.
						//  With event weighting factored in, that number would be variable 
						//  and nTrigger would need to factor in weighting.
						_h1dPhi[centralityIndex]->fill(deltaPhi,1);
					}
				}
			}
		}

		// After all events are processed this is called.
		//  Here I need to normalize my histogram by dividing by
		//  the number of triggers.
		void finalize()
		{
			for (int i = 0; i < N_CENT_TYPES; ++i){
				// calculate the inverse of the number of triggers
				double scaleFactor = 1./(double)nTrigger[i];
				
				// scale by that factor
				_h1dPhi[i]->scaleW(scaleFactor);
				// Histo1Dptr are just pointers to YODA objects, and their documentation
				// is in the YODA docs.
				// e.g.: https://yoda.hepforge.org/doxy/classYODA_1_1Histo1D.html
			}
		}
	};
	
	// The hook for the plugin system
	//  You need to have this, but with your class name, of course.
	DECLARE_RIVET_PLUGIN(STAR_2014_I1234567);
}
