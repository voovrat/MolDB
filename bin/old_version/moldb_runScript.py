import sys
from moldb import Chain,DefList,buildImplementation,parseXparam
from iterate_over_lists import iterate_over_lists
from getopt import gnu_getopt
import os

#scriptText=''
refsList=[];


def print_refs(data):
	global refsList
	params = dict(data);

	prototype = params.pop('moldb_prototype');
	root = params.pop('moldb_root');
	
	visited={}
	m = buildImplementation(params,prototype,visited);
	
	refsList.append( root+'/'+m.getRef() );
	

def writeScript(output,ref_list,xparam_arg,xparam_arg_val):
	scriptText =  "#!/bin/bash\n"	
	if xparam_arg:
		scriptText += '#\n#  XPARAM: '+xparam_arg_val + '\n#\n'

	scriptText += 'FOLDERS="\n'
	
	for ref in ref_list:
		scriptText += ref + '\n';

	scriptText += '"\n'
	scriptText += '\n'
	scriptText += 'N=$(echo $FOLDERS | wc -w )\n'
	scriptText += 'count=0\n'
	scriptText += 'P=$(pwd)\n'
	scriptText += 'for FOLDER in $FOLDERS\n'
	scriptText += 'do\n'
	scriptText += '\techo Molecule $count of $N "(" $((count*100/N)) percent done ")"\n'
	scriptText += '\tcd $FOLDER\n'
	scriptText += '\t. ./PARAMETERS\n'
	scriptText += '\t. ./DEPENDENCIES\n'
	scriptText += '\t./'+script+'\n'
	scriptText += '\tcount=$((count+1))\n'
	scriptText += 'done\n'

	f=open(output,'w');
	f.write(scriptText);
	f.close();

	os.system('chmod +x '+output);
	
	print output


usage="""
Usage: moldb_runScript chain.prototype script [xparam] [-o output] [-x xparam/--xparam=xparams] [-p Npar/--proc Npar]

Generates a script, which will iterate over all the parameters combinations for methodID and runs the script in the corresponding folders

If no iutput filename is specified, the output filename will be generated as:
   PrjectName_ChainName_ScriptName.sh

-x, --xparam=xparams: add/change default parameters of the method.
   xparams in format name1=value1,name2=value2,...
   LIST statements for values are supported (MIND (default) VALUES in LISTS!!!!)

-p, --proc=Npar : number of parallel processes (processors)
  Npar scripts will be generated (+1 main script which runs all others)
"""


if __name__=='__main__':
	if len(sys.argv)<3:
		print usage;
		quit();

	opts,args=gnu_getopt(sys.argv[1:],'o:x:p:',['xparam=','proc='])


#	project_name = args[0];
	chain_fname = args[0];
	chain_name = chain_fname.split('.')[0];

	ch = Chain(chain_fname);
	methodsDict = ch.getMethodsDict();

	project_name = ch.project;
	methodID = ch.runMethod;
#	methodID = args[2];
	prototype = methodsDict[methodID];

	script = args[1];
	script_words=script.split();

	output = project_name + '_' + chain_name + '_'+script_words[0] +'.sh'

	xparam_arg=False
	xparam_arg_val=None
	
	Npar = 1;

	for o in opts:
		k,v = o
		if k=='-o' : output = v;		
		if k=='-x' or k=='--xparam':
			#print 'XPARAM '+v;
			xparam_arg=True
			xparam_arg_val=v
			parseXparam(v,prototype,ch)
			prototype.iterators = prototype.getIterators()
		if k=='-p' or k=='--proc':
			Npar = int(v)

	if len(args)>=3:
		xparam_arg=True
		xparam_arg_val=args[2];
		parseXparam(args[2],prototype,ch)
		prototype.iterators = prototype.getIterators()


	root = '$MOLDB_PROJECTS/'+project_name+'/Methods';

	xparam={'moldb_prototype':prototype,'moldb_root':root };
	iterators = prototype.iterators;


	iterate_over_lists(iterators.keys(),iterators.values(), print_refs, xparam)

	if Npar == 1:
		writeScript(output,refsList,xparam_arg,xparam_arg_val)
	else:
		output_base=output[:-3]
		
		mainText='#!/bin/bash\n';

		N = len(refsList); # / Npar;
	#	beg = 0;

		for i in range(Npar):
			scriptNameBase = output_base + '_p' + str(i)
			scriptName = scriptNameBase +'.sh';
	
			writeScript(scriptName,refsList[i:N:Npar],xparam_arg,xparam_arg_val);

			mainText += './'+scriptName + ' > ' + scriptNameBase +'.out' + ' 2> ' + scriptNameBase + '.err ' +' &\n'
		#	beg += N;
		
		f = open(output,'w');
		f.write(mainText);
		f.close();

	os.system('chmod +x '+output);
	
	



