/*
 * A Root macro to test the functionality of the TreeUtility class
 * 
 */
{
	gSystem->Load("libCintex.so");
 
	Cintex::Enable();
 
	gSystem->Load("libRecoParticleFlowPFClusterTools.so");
	
	TFile f("singleParticle.root");
 
	using namespace std;
	
	using namespace pftools;
	
	TreeUtility tu;
	
	tu.recreateFromRootFile(f);
	
}
