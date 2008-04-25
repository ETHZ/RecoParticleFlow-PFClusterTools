/*
 * A Root macro to test the functionality of the TreeUtility class
 * 
 */
{
	gSystem->Load("libCintex.so");
 
	Cintex::Enable();
 
	gSystem->Load("libRecoParticleFlowPFClusterTools.so");
	
	TFile f("treeUtilityTest.root");
 
	using namespace std;
	
	using namespace minimiser;
	
	TreeUtility tu;
	
	//tu.recreateFromRootFile(f);
	
}
