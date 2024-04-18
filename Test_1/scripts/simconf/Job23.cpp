// Mi headers
#include "/pbs/home/m/mpetro/PROGRAMS/MiModule/include/MiEvent.h" 
#include "/pbs/home/m/mpetro/PROGRAMS/MiModule/include/MiSDVisuHit.h" 
#include "/pbs/home/m/mpetro/PROGRAMS/MiModule/include/MiFilters.h"

#include "TLatex.h"
#include "TVector3.h"


#include <string>
#include <iostream>
#include <sstream>

R__LOAD_LIBRARY(/pbs/home/m/mpetro/PROGRAMS/MiModule/lib/libMiModule.so);

////////////// Function used in script
/////////////////////////////////////////////////////////////

TVector3* get_vertex_vector(MiEvent*  _eve, string _position, int _trackID) // returns the step position of the hit_step when it first enters tracking gas
{
	TVector3* vertexVector;
	if ( _position == "calo" )
	{
		for(int j = 0;j < _eve->getPTD()->getpart(_trackID)->getvertexv()->size();j++)
		{
			if(
				_eve->getPTD()->getpart(_trackID)->getvertex(j)->getpos() == "xcalo" || 
				_eve->getPTD()->getpart(_trackID)->getvertex(j)->getpos() == "calo"  ||
				_eve->getPTD()->getpart(_trackID)->getvertex(j)->getpos() == "gveto" 
			)
			{
				vertexVector = new TVector3(
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getX(),
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getY(),
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getZ()
				);
			}
		}
	}
	else
	{
		for(int j = 0;j < _eve->getPTD()->getpart(_trackID)->getvertexv()->size();j++)
		{
			if(_eve->getPTD()->getpart(_trackID)->getvertex(j)->getpos() == _position)
			{
				vertexVector = new TVector3(
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getX(),
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getY(),
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getZ()
				);
			}
		}
	}

	return vertexVector; //return -1 if the particle never left source foil (happens when we get "fakeItTillYouMakeIt events")
}


float_t get_distance(TVector3* v1, TVector3* v2)
{
	float_t xDifSquared = pow(v1->X() - v2->X(), 2);  // (x1 - x2)^2
	float_t yDifSquared = pow(v1->Y() - v2->Y(), 2);  // (y1 - y2)^2
	float_t zDifSquared = pow(v1->Z() - v2->Z(), 2);  // (z1 - z2)^2

	float_t distance = sqrt( xDifSquared + yDifSquared + zDifSquared );

	return distance; // sqrt( x^2 + y^2 + z^2 )
}

bool is_same_calo_gid(MiGID* cdGID,  MiGID* sdGID )
{
	if(
		cdGID->gettype() ==  sdGID->gettype() &&
		cdGID->getmodule() ==  sdGID->getmodule() &&
		cdGID->getside() ==  sdGID->getside() &&
		cdGID->getwall() ==  sdGID->getwall() &&
		cdGID->getcolumn() ==  sdGID->getcolumn() &&
		cdGID->getrow() ==  sdGID->getrow()
	)
	{
		return true;
	}
	return false;
}

////////////// MAIN BLOCK OF CODE
/////////////////////////////////////////////////////////////
void Job23()
{
////////////// Initialize File names/paths
/////////////////////////////////////////////////////////////
	const char* inFileName                 = "Default.root";													 //FOR TESTING PURPOSES USING ONLY FOLDER 0/
	const char* outPathAngularCorrelation  = "EnePhiDist_Job23.root";
	TFile* 	    outFileAngularCorrelation  = new TFile(outPathAngularCorrelation, "RECREATE");

	int nFiles = 1;

////////////// Initialize variables to be saved
/////////////////////////////////////////////////////////////
	float_t   phi, p1XEscaped, p1YEscaped, p1ZEscaped, p2XEscaped, p2YEscaped, p2ZEscaped;
	float_t   x1Escaped, y1Escaped, z1Escaped, x2Escaped, y2Escaped, z2Escaped;


	TVector3 p1Escaped;
	TVector3 p2Escaped;

	float_t   reconstructedEnergy1, reconstructedEnergy2;
	float_t   simulatedEnergy1, simulatedEnergy2;
	float_t   trackLength1, trackLength2;

////////////// Saving Data
/////////////////////////////////////////////////////////////
	TTree* tree 			= new TTree("tree","tree");

	tree->Branch("phi", &phi, "phi/f");

	tree->Branch("reconstructedEnergy1", &reconstructedEnergy1, "reconstructedEnergy1/f");
	tree->Branch("reconstructedEnergy2", &reconstructedEnergy2, "reconstructedEnergy2/f");

	tree->Branch("simulatedEnergy1", &simulatedEnergy1, "simulatedEnergy1/f");
	tree->Branch("simulatedEnergy2", &simulatedEnergy2, "simulatedEnergy2/f");

	tree->Branch("trackLength1", &trackLength1, "trackLength1/f");
	tree->Branch("trackLength2", &trackLength2, "trackLength2/f");

	tree->Branch("x1Escaped", &x1Escaped, "x1Escaped/f");
	tree->Branch("y1Escaped", &y1Escaped, "y1Escaped/f");
	tree->Branch("z1Escaped", &z1Escaped, "z1Escaped/f");
	tree->Branch("x2Escaped", &x2Escaped, "x2Escaped/f");
	tree->Branch("y2Escaped", &y2Escaped, "y2Escaped/f");
	tree->Branch("z2Escaped", &z2Escaped, "z2Escaped/f");


////////////// Initialize counters
/////////////////////////////////////////////////////////////
	int stepBeforeGas 				= -1;   // Represents the step just before exitting to the tracker gas volume
	int stepBeforeOM 				= -1;

	int nPassed = 0;
	double nPassedFraction = 0.0;

////////////// Loop over files
/////////////////////////////////////////////////////////////
	for( int file = 0; file < nFiles; file ++)
	{
		stringstream ssInPath;
		ssInPath << inFileName;

		if(gSystem->AccessPathName(ssInPath.str().c_str())) // check whether Default.root exists
		{
	    	cout << "Default.root DOESNT EXIST - PATH: " << ssInPath.str().c_str() << endl;
		} 
		else 
		{
			cout << ssInPath.str().c_str() << endl;

	////////////// Initialize reading data
	/////////////////////////////////////////////////////////////
			TFile* 	  		inFile 						= new TFile(ssInPath.str().c_str());
			TTree* 			s 		  					= (TTree*) inFile->Get("Event");
			MiEvent*  		eve 						= new MiEvent();

			s->SetBranchAddress("Eventdata", &eve);
			int nEntries = s->GetEntries();

	/////////////////////////////////////////////////////////////
			for (int e = 0; e < nEntries; ++e) 
			{ 	
				s->GetEntry(e);
					
				cout << " event number = " << e << endl;

				nPassed += 1;
				nPassedFraction += 1.0/double(nEntries);

				reconstructedEnergy1	= eve->getCD()->getcalohit(0)->getE();
				reconstructedEnergy2	= eve->getCD()->getcalohit(1)->getE();

				simulatedEnergy1 = 0;
				simulatedEnergy2 = 0;

				for ( auto & SDCaloHit : *eve->getSD()->getcalohitv() )
				{
					if( is_same_calo_gid(eve->getCD()->getcalohit(0)->getGID(),  SDCaloHit.getGID() ) )
					{
						simulatedEnergy1	+= SDCaloHit.getE(); // there are sometimes multiple SD hits where one will be large and then a few very small (fraction of MeV)
					}
					if( is_same_calo_gid(eve->getCD()->getcalohit(1)->getGID(),  SDCaloHit.getGID() ) )
					{
						simulatedEnergy2	+= SDCaloHit.getE();
					}
				}
				

				TVector3* r1Escaped = get_vertex_vector(eve, "source foil", 0);  	// position vector of the foil vertex
				TVector3* r2Escaped = get_vertex_vector(eve, "source foil", 1);
				TVector3* r1AtOM = get_vertex_vector(eve, "calo", 0);     			// position vector where electron hits OM. Tracklength is calculated as sqrt(r1Escaped^2 + r1AtOM^2)
				TVector3* r2AtOM = get_vertex_vector(eve, "calo", 1);  

				x1Escaped = r1Escaped->X();
				y1Escaped = r1Escaped->Y();
				z1Escaped = r1Escaped->Z();
				x2Escaped = r2Escaped->X();
				y2Escaped = r2Escaped->Y();
				z2Escaped = r2Escaped->Z();

				p1XEscaped = r1AtOM->X() - r1Escaped->X(); // x-coordinate of vector pointing in the direction of electron's travel
				p1YEscaped = r1AtOM->Y() - r1Escaped->Y();
				p1ZEscaped = r1AtOM->Z() - r1Escaped->Z();
				p2XEscaped = r2AtOM->X() - r2Escaped->X();
				p2YEscaped = r2AtOM->Y() - r2Escaped->Y();
				p2ZEscaped = r2AtOM->Z() - r2Escaped->Z();

				p1Escaped = TVector3(p1XEscaped, p1YEscaped, p1ZEscaped);  // This is implemented in MiModule by MP 
				p2Escaped = TVector3(p2XEscaped, p2YEscaped, p2ZEscaped);  // Momentum is returned as TVector3

				phi   = p1Escaped.Angle(p2Escaped)*180/TMath::Pi();

				trackLength1 = get_distance(r1Escaped, r1AtOM);
				trackLength2 = get_distance(r2Escaped, r2AtOM);

				tree->Fill();

			} 	
		}
	}
	// cout << "fakeItTillYouMakeItCounter = " << fakeItTillYouMakeItCounter << endl; // these events are skipped in the simulation 
	cout <<" nPassed = " << nPassed << endl;
	cout <<" nPassedFraction = " << nPassedFraction*100 << "%" << endl;

	outFileAngularCorrelation->Write();
	outFileAngularCorrelation->Close();
}

