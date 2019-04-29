import getopt, sys, subprocess, datetime, gzip, random, os, signal, re


def timenow():
	t = datetime.datetime.today()
	return str(t.year) + "-" + str(t.month).zfill(2) + "-" + str(t.day).zfill(2) + " " + str(t.hour).zfill(2) + ":" + str(t.minute).zfill(2) + ":" + str(t.second).zfill(2)


def FindAFIndex(FormatStrin):
	t = FormatStrin.split(';')
	for x in range(0, len(t)):
		if t[x][:2]=="AF":
			return x
	return -1
	

def FindGTIndex(FormatStrin):
	t = FormatStrin.split(':')
	for x in range(0, len(t)):
		if t[x]=="GT":
			return x
	print ("Format tag GT NOT found in :"+FormatStrin)
	sys.exit()
	

def timediff(t1,t2):
	t = t2-t1
	days = t.days
	secs = t.seconds
	txt = ""
	if days > 0: # at least one day
		txt = txt + str(days) + " days, "
	if secs/3600.0 >= 1: # at least one hour
		minsecs = secs % 3600
		hours = (secs - minsecs)/3600
		txt = txt + str(hours) + " hours, "
		secs = minsecs
	if secs/60.0 >= 1:
		minsecs = secs % 60
		mins = (secs - minsecs)/60
		txt = txt + str(mins) + " minutes, "
		secs = minsecs
	txt = txt + str(secs) + " seconds"
	return txt


def main():
	buildMap = ""
	infile = ""
	outfile = ""
	parsed = 0
	try:
		opts, args = getopt.getopt(sys.argv[1:], "m:v:o:hp", ["mapper=", "vcf=", "out=", "help", "parsed"])
	except getopt.GetoptError, err:
		print str(err)
		sys.exit(2)
	for o,a in opts:
		if o in ["-m","--mapper"]:
			buildMap = a
		elif o in ["-v","--vcf"]:	
			infile = a
		elif o in ["-o","--out"]:	
			outfile = a
		elif o in ["-p","--parsed"]:
			parsed=1
		elif o in ["-h","--help"]:
			usage()
			sys.exit()
		else:
			usage()
			print("Ignoring unknown option: " + o + " " + a)
			sys.exit()

	RanVal = random.randrange(1000,9999)

	if "" in [infile, buildMap]:
		print("\n VcfLiftOver.py\n Parameters in effect:")
		print(" Input VCF file:         [ " + infile + " ]")
		print(" Mapper file:            [ " + buildMap + " ]")
		print("\n [ERROR:] Missing required parameter, exiting.")
		usage()
		sys.exit()
			
	# Start tracking time
	timeTxt = timenow()
	print("Started run at " + timeTxt + "\n\n")
	starttime = datetime.datetime.now()
	
	# Get contig ID from VCF Mapper File metaLines	
	# If Mapper File is VCF format, then parse it to get simplified mapper file
	# if parsed==0:
		# text="zcat "+buildMap+" | head -n 200 | grep -w contig > "+outfile+"."+str(RanVal)+".NEW_CONTIG"
		# subprocess.check_output(text, shell=True, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
		# text=" zcat "+buildMap+" | cut -f 1,2,4,5 | grep -v '#' |  perl -lane 'print \"$F[0]_$F[1]_$F[2]_$F[3]\"' > "+outfile+"."+str(RanVal)+".MAPPER1"
		# os.system(text)
		# text=" zcat "+buildMap+" | cut -f 8 | grep -v '#' | sed '1d' | tr ';' '\\t' | rev | perl -lane 'print \"$F[0]\"' | rev | tr '=' '\\t' | cut -f 2 > "+outfile+"."+str(RanVal)+".MAPPER2"
		# os.system(text)
		# text=" paste  "+outfile+"."+str(RanVal)+".MAPPER1 " +outfile+"."+str(RanVal)+".MAPPER2 > "+outfile+"."+str(RanVal)+".MAPPER"
		# os.system(text)
		# buildMap = outfile+"."+str(RanVal)+".MAPPER"
	# else:
		# print "\nUsing user-input parsed Map file ...\n"
	dict = {}

	# Read and store list of variants from mapper file
	with open(buildMap, "r") as f:
		for line in f:
			values = line.rstrip("\n").split('\t')
			dict[values[0]] = values[1]

	#with gzip.open(outfile+".vcf.gz","w") as wm: 
	
	Match=0
	Swap=0
	Miss=0
	Start=0
	with gzip.open(infile,"r") as f, open(outfile+"."+str(RanVal)+".tempvcf","w") as wm, open(outfile+".missVariants","w") as miss:
		for line in f:
			temp = line.rstrip("\n").split()
			IsIndic = 0
			if line[0]=='#':
				wm.write(line)
				if parsed==0:
					if Start==0:
						with open("/net/fantasia/home/sayantan/DATABASE/BUILD38_VARIANT_MAPPER/build38.contig.header","r") as con:
							for lineCont in con:
								wm.write(lineCont)
				Start=1		 
			else:
				VarName = temp[0]+"_"+temp[1]+"_"+temp[3]+"_"+temp[4]
				SwapVarName = temp[0]+"_"+temp[1]+"_"+temp[4]+"_"+temp[3]
				if VarName in dict:
					IsIndic=1
				elif SwapVarName in dict:
					Swap=Swap+1
					IsIndic=2
				
				if IsIndic>0:
					Match=Match+1
					GTIndex = FindGTIndex(temp[8])
					AFIndex = FindAFIndex(temp[7])
					if AFIndex>-1:
						AFVal = ((temp[7].split(";"))[AFIndex]).split("=")[1]
						if IsIndic==1:
							temp[7]="AF="+AFVal
						else:
							temp[7]="AF="+str(1-float(AFVal))
							#print SwapVarName+"\t"+AFVal+"\t"+temp[7]
					else:
						temp[7]="."
				
					if IsIndic==1:
						tempNew=dict[VarName].split('_')
					else:
						tempNew=dict[SwapVarName].split('_')

					temp[0]=tempNew[0]
					temp[1]=tempNew[1]
					if IsIndic==2:
						temp[3]=tempNew[2]
						temp[4]=tempNew[3]
					
						
					temp[8]="GT"
					for x in range(9, len(temp)):
						GTVal = (temp[x].split(":"))[GTIndex]
						GTVal=GTVal.replace('/','|')
						GT = GTVal.split("|")
						
						#print VarName
						
						if IsIndic==2:
							if GT[0]!='.':
								GT[0]=str(1-int(GT[0]))
							if GT[1]!='.':
								GT[1]=str(1-int(GT[1]))
						temp[x]=("|".join(GT))
					wm.write("\t".join(temp))
					wm.write("\n")
				else:	
					Miss=Miss+1
					miss.write(VarName+"\n")

	os.system("/usr/cluster/bin/vcf-sort "+outfile+"."+str(RanVal)+".tempvcf |  /usr/cluster/bin/bgzip -c > "+outfile+".vcf.gz")
	os.system("/usr/cluster/bin/tabix -p vcf "+outfile+".vcf.gz")
	
	# Put at the end of run
	endtime = datetime.datetime.now()
	timeTxt = timenow()
	print("\n "+ str(Match) + " variants lifted over (" +str(Swap) + " variants had REF/ALT swapped)")
	print("\n "+ str(Miss) + " variants skipped (listed in "+ outfile+".missVariants)")
	os.system("rm "+outfile+"."+str(RanVal)+".*")
	print("\nFinished run at " + timeTxt)
	print("Total Runtime: " + timediff(starttime, endtime))
	

def usage():
	print(''' \n VcfLiftOver.py -- Lift Over VCF file using input mapper
	
 Usage: python VcfLiftOver.py -m [Mapper.txt] -v [Input.vcf.gz] -o [Ouptut.vcf.gz] -p 

 Required parameters:
  -v [Input.vcf.gz] : Input VCF file
  -m [Mapper.txt]   : VCF file containing build lift over information
  -o [Ouptut.vcf.gz]: Output VCF file
  -p                : Mapper File is already parsed as follows:  File
                      containing current build variants in first column 
                      and new build variants in the second column. Format for
                      variants should be [CHR_POS_REF_ALT]
  Written by Sayatan Das, sayantan@umich.edu
  Last edited 2017-09-12
''')


def version():
	print('''Version notes:

v0.1-2017-09-12: First release. Does not handle multiallelics.
''')


if __name__ == "__main__":
	main()
