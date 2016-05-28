#!/usr/bin/env julia0.3.11

########
# script analysing the information available in the multifasta files
# This version is able to distingush between the whole dataset and the ingroup
# but is HORRIBLY SLOW!!!
########

## I) decalre data type and functions
# data type not strictly needed, but maybe useful if we want to extend later
type Chunk
	# individual name: a string
	race::String
	# individual name: a string
	nom::String
	# haplotype as integer
	haplo::Int
	# sequence: a list of characters
	data::Vector{Char}
	
end

# add race name
function add_race(chunk, spline)
	race = split(spline[3], '_')[1]
	chunk.race = race
end

# add sample name
function add_nom(chunk, spline)
	nom = spline[3]
	chunk.nom = nom
end

# add haplotype number
function add_haplo(chunk, spline)
	hap = parseint(spline[4])
	chunk.haplo = hap
end

# add a line of the sequence
function add_data(chunk, line)
	# collect converts string into array of char
	append!(chunk.data, collect(line))
end

# read the entire file
function read_fasta(input)
	# 1 chunk = 1 individual
	chunks = Chunk[]

	for line in eachline(input)
#	for line in input
		# strip removes whitespace from beginning and end
		line = strip(line)

		# ignore empty lines
		if isempty(line)
			continue
		end

		# header
		if line[1] == '>'
			# create a chunk from the line
			spline = split(line, [' ', '#', ':'], 0, false)
			push!(chunks, Chunk("", "", -1, Char[]))
			add_race(chunks[end], spline)
			add_nom(chunks[end], spline)
			add_haplo(chunks[end], spline)
		# data
		else
			add_data(chunks[end], line)
		end
	end
	
	# last expression -> return value of function
	chunks
end

# check information at each site and update corresponding vectors
function check_info(i, ref, info, hasN, infoN, hasIns, infoNIns, infoIns)
	if ref[i] == 'A' || ref[i] == 'T' ||
		ref[i] == 'G' || ref[i] == 'C'
	
		info[i] = 1
	
		# 1.4.0.1) informative with N
		if hasN[i] == 1
			infoN[i] = 1
		
			# 1.4.0.1.1) informative with N and indels
			if hasIns[i] == 1
				infoNIns[i] = 1
			end
		end
	
		# 1.4.0.2) informative with indels
		if hasIns[i] == 1
			infoIns[i] = 1
		end
	end
	info, infoN, infoIns, infoNIns
end


# check type of polymorphism
function check_polym(i, hasN, hasIns, polCle, polN, polNIns, polIns)
	# 1.4.2.1) count the number of clean polym
	if hasN[i] == 0 && hasIns[i] == 0
		polCle[i] = 1
	else
		# 1.4.2.2) count number of polymorphic sites with N
		if hasN[i] == 1
			polN[i] = 1
		
			# 1.4.2.2.1) nber of polymorphic sites with N and indels
			if hasIns[i] == 1
				polNIns[i] = 1
			end
		end
	
		# 1.4.2.3) count number of polymorphic sites with indels
		if hasIns[i] == 1
			polIns[i] = 1
		end
	end
	polCle, polN, polNIns, polIns
end


# calculate stats
function stats(chunks::Vector{Chunk}, vecfile, fasta, outgroups)

	outg = strip(outgroups)
	ref = []
	refIng = []

	# whole dataset
		# list of 0s; we'll need int later anyway
	hasN = fill(0, length(chunks[1].data))
	countN = 0
	
	hasIns = fill(0, length(chunks[1].data))
	countIns = 0
	
	info = fill(0, length(chunks[1].data))
	infoClean = fill(0, length(chunks[1].data))
	infoN = fill(0, length(chunks[1].data))
	infoIns = fill(0, length(chunks[1].data))
	infoNIns = fill(0, length(chunks[1].data))
	
	polym = fill(0, length(chunks[1].data))
	polCle = fill(0, length(chunks[1].data))
	polN = fill(0, length(chunks[1].data))
	polIns = fill(0, length(chunks[1].data))
	polNIns = fill(0, length(chunks[1].data))

	# Ingroup
		# list of 0s; we'll need int later anyway
	hasNIng = fill(0, length(chunks[1].data))
	countNIng = 0
	
	hasInsIng = fill(0, length(chunks[1].data))
	countInsIng = 0
	
	infoIng = fill(0, length(chunks[1].data))
	infoCleanIng = fill(0, length(chunks[1].data))
	infoNIng = fill(0, length(chunks[1].data))
	infoInsIng = fill(0, length(chunks[1].data))
	infoNInsIng = fill(0, length(chunks[1].data))
	
	polymIng = fill(0, length(chunks[1].data))
	polCleIng = fill(0, length(chunks[1].data))
	polNIng = fill(0, length(chunks[1].data))
	polInsIng = fill(0, length(chunks[1].data))
	polNInsIng = fill(0, length(chunks[1].data))

	incIng = 0	# switch ref sequence ingoup
	for c in 1:length(chunks)
		chunk = chunks[c]
		rac = chunk.race

		# make sure all chunks are the same length
		@assert	length(chunk.data) == length(hasN)
		
		# 0) set ref sequences
			# whole dataset
		if c == 1
			ref = chunk.data
		end
		
			# ingroup
		if incIng == 0 && (!ismatch(Regex(rac), outg))
			refIng = chunk.data
			incIng += 1
		end
		
		# start incrementation individual informative sites
		indi_info = 0
		
		# 1) check all sites
		for i in 1:length(chunk.data)

			# 1.1) count number of Ns
			if chunk.data[i] == 'N'
				# whole dataset
				countN += 1
				hasN[i] = 1
			
				# ingroup
				if !ismatch(Regex(rac), outg)
					countNIng += 1
					hasNIng[i] = 1
				end
			

			# 1.2) count number of indel positions
			elseif chunk.data[i] == '-'
				# whole dataset
				countIns += 1
				hasIns[i] = 1
				
				# ingroup
				if !ismatch(Regex(rac), outg)
					countInsIng += 1
					hasInsIng[i] = 1
				end

			# 1.3) record informative sites
			elseif chunk.data[i] == 'A' || chunk.data[i] == 'T' || 
				chunk.data[i] == 'G' || chunk.data[i] == 'C'
				
				# 1.3.1) count number of informative sites for individual i
				indi_info += 1

				# 1.3.2) count number of polymorphic sites
					# 1.3.2.1) update ref sequences
						# whole dataset
				if ref[i] == 'N' || ref[i] == '-'
					ref[i] = chunk.data[i]
				end
				
						# ingroup
				if (!ismatch(Regex(rac), outg)) && (refIng[i] == 'N' || refIng[i] == '-')
					refIng[i] = chunk.data[i]
				end

					# 1.3.2.2) check polymorphism
				if polymIng[i] == 0 && chunk.data[i] != ref[i]
					polym[i] = 1	# whole dataset
					if (!ismatch(Regex(rac), outg)) && chunk.data[i] != refIng[i]
						polymIng[i] = 1	# ingroup
					end
				end
				
			else
				error("Invalid symbol" * chunk.data[i] * " in alignment: file "
					 * fasta * ", haplotype " * c * " (valid symbols are 'A', 'C', 'G', 'T', 'N' and '-')")
			end

			#### print stats for each individual ####
			if i == length(chunk.data)
				println(vecfile, fasta, "\t", chunk.race, "\t", chunk.nom, "\t",  chunk.haplo,"\t",  indi_info) 
				indi_info = 0
			end
			#### print stats for each individual ####


			# 1.4) if last chunk, perform few more tests
			if c == length(chunks)

				# 1.4.1) check the information of the site i
					# 1.4.1.1) clean?
				if hasN[i] == 0 && hasIns[i] == 0	# automatically info
					infoClean[i] = 1	# whole
					infoCleanIng[i] = 1	# ingroup
				else
					# 1.4.1.2) check info whole dataset
				info, infoN, infoIns, infoNIns = check_info(i, ref, info, hasN, infoN, hasIns, infoNIns, infoIns)

					# 1.4.1.3) check if ingroup clean
					if hasNIng[i] == 0 && hasInsIng[i] == 0
						infoCleanIng[i] = 1
					else
						infoIng, infoNIng, infoInsIng, infoNInsIng = check_info(i, refIng, infoIng, hasNIng, infoNIng, hasInsIng, infoNInsIng, infoInsIng)
					end
				end

				# 1.4.2) for polymorphic sites
				if polymIng[i] == 1
					# whole dataset
					polCle, polN, polNIns, polIns = check_polym(i, hasN, hasIns, polCle, polN, polNIns, polIns)

					# ingroup
					polCleIng, polNIng, polNInsIng, polInsIng = check_polym(i, hasNIng, hasInsIng, polCleIng, polNIng, polNInsIng, polInsIng)
				else
					if  polym[i] == 1
						# whole dataset
					polCle, polN, polNIns, polIns = check_polym(i, hasN, hasIns, polCle, polN, polNIns, polIns)
					end
				end

			end	# end of 'if c == length(chunks)'
		end	# end of 'for i in 1:length(chunk.data)'
	end	# end of 'for c in 1:length(chunks)'

	# functions can return several values
	propN = countN/(length(hasN)*length(chunks))
	propIns = countIns/(length(hasIns)*length(chunks))
	propNIng = countNIng/(length(hasNIng)*length(chunks))
	propInsIng = countInsIng/(length(hasInsIng)*length(chunks))
	
	hasN, propN, hasIns, propIns, infoClean, polym, polCle, polN, polIns, polNIns, info, infoN, infoIns, infoNIns, hasNIng, propNIng, hasInsIng, propInsIng, infoCleanIng, polymIng, polCleIng, polNIng, polInsIng, polNInsIng, infoIng, infoNIng, infoInsIng, infoNInsIng
end

function print_help()
	println("--------------")
	println("usage:")
	println("fastaNstats.jl [-o OUTRGOUP1[,OUTGROUP2,...]] -s <INDIV_STAT_FILE> -f <FASTA_FILES>...\n")
	println("WARNING: <INDIV_STAT_FILE> cannot cannot end with any of the usual fasta file extensions (see 'http://en.wikipedia.org/wiki/FASTA_format#File_extension')")
	println("--------------")
end

function get_arg(args, i, fun = x->x)
	if length(args) < i+1
		error("expected argument after $(args[i])")
	end
	i += 1
	fun(args[i]), i
end

## II) script
if length(ARGS) < 2 || ismatch(r"\.fasta$|\.fas$|\.fa$|\.seq$|\.fsa$|\.fna$|\.ffn$|\.faa$|\.frn$", ARGS[1])
	print_help()
	exit()
end

type Args
	outgroups
	indivstat
	j
end

args = Args("", "", -1)

i = 1
while i ≤ length(ARGS)
	arg = ARGS[i]

	if arg == "-o"
		args.outgroups, i= get_arg(ARGS, i)
		println("# Outgroups are $(args.outgroups)")
	elseif arg == "-s"
		args.indivstat, i = get_arg(ARGS, i)
		println("# File to print individual stats is $(args.indivstat)")
	elseif arg == "-f"
		args.j = i + 1
		println("# The first fasta file to analyse is $(ARGS[args.j])")
	end
	i += 1
end

# open INDIV_STAT_FILE file
vecfile = open(args.indivstat, "w")
println(vecfile, "file\trace\tclone\thaplotype\tNbInfoSites")

# print header of stdout
println("fasta\tlg\tNpos\tNposIng\tPcN\tPcNIng\tNinfo\tNinfoIng\tNinfoClean\tNinfoCleanIng\tNinfoN\tNinfoNIng\tNinfoIns\tNinfoInsIng\tNinfoNIns\tNinfoNInsIng\tNpolym\tNpolymIng\tNpolCle\tNpolCleIng\tNpolN\tNpolNIng\tNpolIns\tNpolInsIng\tNpolNIns\tNpolNInsIng\tpropN\tpropNIng\tInspos\tInsposIng\tpropInsPos\tpropInsPosIng\tpropIns\tpropInsIng\tfullColN")

for f in args.j:length(ARGS)

	fasta = ARGS[f]

	# apply read_fasta to file, then close it again
	chunks = open(read_fasta, fasta)

	hasN, propN, hasIns, propIns, infoClean, polym, polCle, polN, polIns, polNIns, info, infoN, infoIns, infoNIns, hasNIng, propNIng, hasInsIng, propInsIng, infoCleanIng, polymIng, polCleIng, polNIng, polInsIng, polNInsIng, infoIng, infoNIng, infoInsIng, infoNInsIng = stats(chunks, vecfile, fasta, args.outgroups)

	lg = length(hasN)
	
	# whole dataset
	Ninfo = sum(info)
	NinfoClean = sum(infoClean)
	NinfoN = sum(infoN)
	NinfoIns = sum(infoIns)
	NinfoNIns = sum(infoNIns)
	
	Npolym = sum(polym)
	NpolCle = sum(polCle)
	NpolN = sum(polN)
	NpolIns = sum(polIns)
	NpolNIns = sum(polNIns)
	
	Npos = sum(hasN)
	PcN = Npos/lg
	fullColN = propN == PcN

	Inspos = sum(hasIns)
	propInsPos = Inspos/lg
	
	# ingroup
	NinfoIng = sum(infoIng)
	NinfoCleanIng = sum(infoCleanIng)
	NinfoNIng = sum(infoNIng)
	NinfoInsIng = sum(infoInsIng)
	NinfoNInsIng = sum(infoNInsIng)
	
	NpolymIng = sum(polymIng)
	NpolCleIng = sum(polCleIng)
	NpolNIng = sum(polNIng)
	NpolInsIng = sum(polInsIng)
	NpolNInsIng = sum(polNInsIng)
	
	NposIng = sum(hasNIng)
	PcNIng = NposIng/lg

	InsposIng = sum(hasInsIng)
	propInsPosIng = InsposIng/lg

	println(fasta, "\t", lg, "\t",
			Npos, "\t", NposIng, "\t",
			PcN, "\t", PcNIng, "\t",
			Ninfo, "\t", NinfoIng, "\t",
			NinfoClean, "\t", NinfoCleanIng, "\t",
			NinfoN, "\t", NinfoNIng, "\t",
			NinfoIns, "\t", NinfoInsIng, "\t",
			NinfoNIns, "\t", NinfoNInsIng, "\t",
			Npolym, "\t", NpolymIng, "\t",
			NpolCle, "\t", NpolCleIng, "\t", 
			NpolN, "\t", NpolNIng, "\t",
			NpolIns, "\t", NpolInsIng, "\t", 
			NpolNIns, "\t", NpolNInsIng, "\t",
			propN, "\t", propNIng, "\t",
			Inspos, "\t", InsposIng, "\t",
			propInsPos, "\t", propInsPosIng, "\t",
			propIns, "\t", propInsIng, "\t",
			fullColN)
end

# close outfile
close(vecfile)

