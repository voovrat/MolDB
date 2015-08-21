from moldb import SwitchChain
import sys

usage="""
moldb_generateSwitches  File.switch

creates several chain files (without SWITCHes) using one file with SWITCHes
The names of files are determined by the names of CASEs
"""

if __name__ == '__main__':
	if len(sys.argv)<2:
		print usage
		quit();

	sw = SwitchChain(sys.argv[1]);
	chains,chain_names = sw.generateChains();

	d = sw.getMethodsDict();

	MethodName = d[sw.runMethod].name;

	for chid in range(len(chains)):
		ch = chains[chid];
		chname = chain_names[chid];
	
		fname = MethodName + chname + '.chain';

		print fname

		f = open(fname,'w');
		f.write(str(ch));
		f.close();


