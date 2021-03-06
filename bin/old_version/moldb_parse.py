from moldb import skipEmptyLines,raiseException,DefList,str2val

#
# Parses SWITCH statement.
# Returns the list of all different combinations of parameters
#
# EXAMPLE:  SWITCH
#			2 CASES
#				CASE 1
#				2 RECORDS
#					a = LIST 2 5 4 (default)
#					b = FILE a.txt
#				ENDCASE
#
#				CASE 2
#				1 RECORD
#					c = LIST 2 (default) 8 6
#				ENDCASE
#			ENDSWITCH
#
#		will return [ {a: LIST 2 5 4, b: LIST a.txt1rec atxt2rec...},{c:LIST 2 8 6}]
#		here LIST means DefList (see moldb.py)
#
def parseSwitch(lines,currline):
	
	currline = skipEmptyLines(lines,currline);
	
	# SWITCH
	words = lines[currline].split();
	if words[0]!='SWITCH':
		raiseException(currline,'SWITCH expected');

	# n CASES
	currline = skipEmptyLines(lines,currline+1);
	words = lines[currline].split();
	if len(words)<2 or words[1]!='CASES':
		raiseException(currline,'<number> CASES expected');

	try:
		Ncase = int(words[0]);
	except ValueError:
		raiseException(currline,'<number> CASES is expected (but not valid number '+words[0] +' given)');


	params_variety = [];
	param_names=[];

	for i in range(Ncase):
		currline = skipEmptyLines(lines,currline+1);
		currline,sub_param_variety,sub_param_names = parseCase(lines,currline);

		#print 'sub_param_variety:'+str(sub_param_variety)
		# sub_param_variety is a List(!) of parameters, ( because switches can be nested !!! )
		params_variety += sub_param_variety;
		param_names += sub_param_names;

	currline = skipEmptyLines(lines,currline+1);

	words = lines[currline].split();
	if words[0]!='ENDSWITCH':
		raiseException(currline,'ENDSWITCH expected');

	return currline,params_variety,param_names

#
# Parses statement 
#  CASE  <name>
#     PArameters
# ENDCASE
# 
# returns a param variety
#
def parseCase(lines,currline):
	currline = skipEmptyLines(lines,currline);
	words = lines[currline].split();

	if len(words)<2 or words[0]!='CASE':
		raiseException(currline,'CASE <Name> is expected');

	case_name = words[1];

	currline=skipEmptyLines(lines,currline+1);

	currline,param_variety,param_names = parseParameters(lines,currline)

	new_names = [];

	for name in param_names:
		new_names.append(case_name+name);

	param_names = new_names;

	currline=skipEmptyLines(lines,currline+1);

	words=lines[currline].split();

	if words[0]!='ENDCASE':
		raiseException(currline,'ENDCASE expected');

	return currline,param_variety,param_names


# Parses the following structure:
#
#  N RECORDS
#
#	name1=value1
#	name2=value2
#
#	SWITCH
#		...
#	ENDSWITCH
#
#	name3=value3
#
#
# NOTE! SWITCH is counted as one record! 
#
#  returns the currline,param_varitely (list of parameters dictionaries for different SWITCH options!)
#
def parseParameters(lines,currline):

	currline = skipEmptyLines(lines,currline);
	
	words=lines[currline].split();

	if len(words)<2 or words[1]!='RECORD' and words[1]!='RECORDS':
		raiseException(currline,'<number> RECORDS expected');

	NRec = int(words[0]);

	param_variety = [{}];
	param_names = [''];

	for i in range(NRec):
		currline = skipEmptyLines(lines,currline+1);

		tokens = lines[currline].split();
		if tokens[0] == 'SWITCH':
			currline,sub_param_variety,sub_param_names = parseSwitch(lines,currline);

			new_variety = [];

			for oldprm in param_variety:
				for newprm in sub_param_variety:
					joinprm = dict(oldprm);
					joinprm.update(newprm);
					new_variety.append(joinprm);
				
					
			param_variety = new_variety;

			new_names = [];
			for oldname in param_names:
				for newname in sub_param_names:
					new_names.append(oldname+':' + newname);

			param_names = new_names;
			
		else:
			words = lines[currline].split('=');

			param = words[0].strip().lower();
		
			values = words[1].split();
			keyword = values[0];

			if keyword=='LIST':
				value = DefList(values[1:]);
			elif keyword=='FILE':
				f=open(values[1]);
				file_values = [ str2val(s.strip()) for s in f.read().splitlines()];
				f.close();

				value = DefList(file_values);
			else:
				value = str2val(values[0]);

			for prm in param_variety:
				prm[param]=value;

	return currline,param_variety,param_names;

if __name__ == '__main__':
	f=open('simple.test');
	lines=f.readlines();
	f.close();

	currline,param_variety,param_names = parseParameters(lines,0);

	print 'RESULT:'+str(param_variety)
	print 'NAMES:'+str(param_names)
