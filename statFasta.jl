#!/usr/bin/env julia0.3.11

########
# script analysing the information available in the multifasta files
# Optimized version of statFasta_2.0.jl (compute stats for the ingroup and/or the whole dataset)
# almost as quick as statFasta_1.0.jl that computes stats only for the whole dataset
########

## I) declare data type and functions
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

# type of the object recording all the stats
	# WARNING: think to modify the 'set_stat' if this type is modified!
type Stats
	hasN::Vector{Int}
	countN::Int
	
	hasIns::Vector{Int}
	countIns::Int
	
	called::Vector{Int}         # sites positively called (with at least one nucleotide)
	calledClean::Vector{Int}    # sites with nucleotides only (no N, no indel)
	calledN::Vector{Int}        # called sites with at least one N
	calledIns::Vector{Int}      # idem with indels
	calledNIns::Vector{Int}     # called sites with at least one N and one indel
	
	polym::Vector{Int}          # polymorphic sites
	polCle::Vector{Int}         # polymorphic clean (no N, no indel)
	polN::Vector{Int}           # polymorphic site with at least one N
	polIns::Vector{Int}         # polymorphic site with at least one indel
	polNIns::Vector{Int}        # polymorphic site with at least one N and one indel
	
	propN::Float64              # proportion of sites with N in the alignement
	propIns::Float64            # proportion of sites with N and indel in the alignement
end


# set stats object
function set_stat(lg)
# where lg is the length of the sequence
	stats = Stats(
		fill(0, lg), 0, fill(0, lg), 0,
		fill(0, lg), fill(0, lg), fill(0, lg), 
		fill(0, lg), fill(0, lg), fill(0, lg), 
		fill(0, lg), fill(0, lg), fill(0, lg), 
		fill(0, lg),
		-1, -1
	)
	stats
end

# add race name
function add_race(chunk, spline)
	race = split(spline[3], '_')[1]
	if race == ""
		error("Incorrect race name in $(spline)")
	else
		chunk.race = race
	end
	chunk.race
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

# check call at each site
function check_called(i, ref, stats)
	# 1) check if site is has been positively called
	if ref[i] == 'A' || ref[i] == 'T' ||
		ref[i] == 'G' || ref[i] == 'C'
	
		stats.called[i] = 1
	
		# 1.1) positively called with N
		if stats.hasN[i] == 1
			stats.calledN[i] = 1
		
			# 1.2) positively called with N and indels
			if stats.hasIns[i] == 1
				stats.calledNIns[i] = 1
			end
		end
	
		# 2) positively called with indels
		if stats.hasIns[i] == 1
			stats.calledIns[i] = 1
		end
	end
	stats
end

# check type of polymorphism
function check_polym(i, stats)
	# 1) count the number of clean polym
	if stats.hasN[i] == 0 && stats.hasIns[i] == 0
		stats.polCle[i] = 1
	else
		# 2.1) count number of polymorphic sites with N
		if stats.hasN[i] == 1
			stats.polN[i] = 1
		
			# 2.2) nber of polymorphic sites with N and indels
			if stats.hasIns[i] == 1
				stats.polNIns[i] = 1
			end
		end
	
		# 2.3) count number of polymorphic sites with indels
		if stats.hasIns[i] == 1
			stats.polIns[i] = 1
		end
	end
	stats
end

# function describe site calling/polymorphism
function descSite (i, ref, whole, refIng, ingrp)
	# 1) describe call of site i
		# 1.1) clean? If so, automatically called and calledClean
	if whole.hasN[i] == 0 && whole.hasIns[i] == 0
        whole.called[i] = 1
		ingrp.called[i] = 1
		whole.calledClean[i] = 1
		ingrp.calledClean[i] = 1
	else
		# 1.2) check if called whole dataset
		whole = check_called(i, ref, whole)

		# 1.3) check if ingroup clean
		if ingrp.hasN[i] == 0 && ingrp.hasIns[i] == 0
			ingrp.calledClean[i] = 1
		else
			ingrp = check_called(i, refIng, ingrp)
		end
	end
	
	# 2) define polymorphism at site i
	if ingrp.polym[i] == 1
		whole = check_polym(i, whole)
		ingrp = check_polym(i, ingrp)
	else
		if  whole.polym[i] == 1
			# whole dataset
		whole = check_polym(i, whole)
		end
	end
	whole, ingrp
end

# calculate stats
function compStats(chunks::Vector{Chunk}, vecfile, fasta, outgroups)

	outg = split(strip(outgroups), ",")
	
	ref = []
	refIng = []

	# set constructors of statistics
		# whole dataset
	whole = set_stat(length(chunks[1].data))
		# Ingroup
	ingrp = set_stat(length(chunks[1].data))

	incIng = 0	# switch ref sequence ingoup
	for c in 1:length(chunks)
		chunk = chunks[c]
		rac = chunk.race
		
		id_outg = find (outg .== rac)
		is_outg = length(id_outg) == 0 ? false : true

		# make sure all chunks are the same length
		@assert	length(chunk.data) == length(whole.hasN)
		
		### 0) set ref sequences
			# whole dataset
		if c == 1
			ref = chunk.data
		end
		
			# ingroup
		if incIng == 0 && !is_outg
			refIng = chunk.data
			incIng += 1
		end
		
		# start incrementation individual positively called sites
		indi_called = 0
		
		### 1) check all sites
		# 1.1) if individual is from the outgroup, update for the whole dataset only
		if is_outg || outg == ""
			# loop over all positions
			for i in 1:length(chunk.data)
				# 1.1.1) count number of Ns
				if chunk.data[i] == 'N'
					whole.countN += 1
					whole.hasN[i] = 1
					
				# 1.1.2) count number of indel positions
				elseif chunk.data[i] == '-'
					whole.countIns += 1
					whole.hasIns[i] = 1
					
				# 1.1.3) record positively called sites
				elseif chunk.data[i] == 'A' || chunk.data[i] == 'T' || 
				chunk.data[i] == 'G' || chunk.data[i] == 'C'
				
					# 1.1.3.1) count number of positively called sites for individual i
					indi_called += 1
					
					# 1.1.3.2) count number of polymorphic sites
						# 1.1.3.2.1) update ref sequences
					if ref[i] == 'N' || ref[i] == '-'
						ref[i] = chunk.data[i]
					end
					
						# 1.1.3.2.2) check polymorphism
					if chunk.data[i] != ref[i]
						whole.polym[i] = 1	# whole dataset
					end
				else
					error("Invalid symbol" * chunk.data[i] * " in alignment: file "
						 * basename(fasta) * ", haplotype " * c * 
						 " (valid symbols are 'A', 'C', 'G', 'T', 'N' and '-')")
				end
				
				# 1.1.4) if last chunk, describe information of the site i
				if c == length(chunks)
					whole, ingrp = descSite (i, ref, whole, refIng, ingrp)
				end
			
			end # end of loop over all positions
			
			# 1.1.5) print per individual stats
			println(vecfile, basename(fasta), "\t", chunk.race, "\t", chunk.nom,
					 "\t",  chunk.haplo,"\t",  indi_called) 


		# 1.2) if individual is from the ingroup, update for the ingrp.dataset plus the ingroup, independently
		else
			# loop over all positions
			for i in 1:length(chunk.data)
				# 1.2.1) count number of Ns
				if chunk.data[i] == 'N'
					whole.countN += 1
					whole.hasN[i] = 1
					ingrp.countN += 1
					ingrp.hasN[i] = 1
					
				# 1.2.2) count number of indel positions
				elseif chunk.data[i] == '-'
					whole.countIns += 1
					whole.hasIns[i] = 1
					ingrp.countIns += 1
					ingrp.hasIns[i] = 1
					
				# 1.2.3) record positively called sites
				elseif chunk.data[i] == 'A' || chunk.data[i] == 'T' || 
				chunk.data[i] == 'G' || chunk.data[i] == 'C'
				
					# 1.2.3.1) count number of positively called sites for individual i
					indi_called += 1
					
					# 1.2.3.2) count number of polymorphic sites
						# 1.2.3.2.1) update ref sequences
					if ref[i] == 'N' || ref[i] == '-'
						ref[i] = chunk.data[i]
					end
					if refIng[i] == 'N' || refIng[i] == '-'
						refIng[i] = chunk.data[i]
					end
					
						# 1.2.3.2.2) check polymorphism
					if chunk.data[i] != ref[i]
						whole.polym[i] = 1	# whole dataset
					end
					if chunk.data[i] != refIng[i]
						ingrp.polym[i] = 1	# ingrp.dataset
					end
					
				else
					error("Invalid symbol" * chunk.data[i] * " in alignment: file "
						 * basename(fasta) * ", haplotype " * c * " (valid symbols are 'A', 'C', 'G', 'T', 'N' and '-')")
				end
			
				# 1.2.4) if last chunk, describe positively called of the site i
				if c == length(chunks)
					whole, ingrp = descSite (i, ref, whole, refIng, ingrp)
				end
			end # end of loop over all positions
		
			# 1.2.5) print per individual stats
			println(vecfile, basename(fasta), "\t", chunk.race, "\t", chunk.nom,
					 "\t",  chunk.haplo,"\t",  indi_called) 

		end # end of analysis of sequence block (outgroup/ingroup) 
		
	end	# end of 'for c in 1:length(chunks)'


	# 2) prepare data for returning
	whole.propN = whole.countN/(length(whole.hasN)*length(chunks))
	whole.propIns = whole.countIns/(length(whole.hasIns)*length(chunks))
	ingrp.propN = ingrp.countN/(length(ingrp.hasN)*length(chunks))
	ingrp.propIns = ingrp.countIns/(length(ingrp.hasIns)*length(chunks))
	
		# functions can return several values
	whole, ingrp
end

function lastStats(whole, ingrp)
	lg = length(whole.hasN)

	# whole dataset
	Ncalled = sum(whole.called)
	NcalledClean = sum(whole.calledClean)
	NcalledN = sum(whole.calledN)
	NcalledIns = sum(whole.calledIns)
	NcalledNIns = sum(whole.calledNIns)
	
	Npolym = sum(whole.polym)
	NpolCle = sum(whole.polCle)
	NpolN = sum(whole.polN)
	NpolIns = sum(whole.polIns)
	NpolNIns = sum(whole.polNIns)
	
	Npos = sum(whole.hasN)
	PcN = Npos/lg
	fullColN = whole.propN == PcN

	Inspos = sum(whole.hasIns)
	propInsPos = Inspos/lg
	propN = whole.propN
	propIns = whole.propIns
	
	# ingroup
	NcalledIng = sum(ingrp.called)
	NcalledCleanIng = sum(ingrp.calledClean)
	NcalledNIng = sum(ingrp.calledN)
	NcalledInsIng = sum(ingrp.calledIns)
	NcalledNInsIng = sum(ingrp.calledNIns)
	
	NpolymIng = sum(ingrp.polym)
	NpolCleIng = sum(ingrp.polCle)
	NpolNIng = sum(ingrp.polN)
	NpolInsIng = sum(ingrp.polIns)
	NpolNInsIng = sum(ingrp.polNIns)
	
	NposIng = sum(ingrp.hasN)
	PcNIng = NposIng/lg

	InsposIng = sum(ingrp.hasIns)
	propInsPosIng = InsposIng/lg
	propNIng = ingrp.propN
	propInsIng = ingrp.propIns
	
	lg, Ncalled, NcalledClean, NcalledN, NcalledIns, NcalledNIns, Npolym, NpolCle, NpolN, NpolIns, NpolNIns, Npos, PcN, fullColN, Inspos, propInsPos, propN, propIns, NcalledIng, NcalledCleanIng, NcalledNIng, NcalledInsIng, NcalledNInsIng, NpolymIng, NpolCleIng, NpolNIng, NpolInsIng, NpolNInsIng, NposIng, PcNIng, InsposIng, propInsPosIng, propNIng, propInsIng
end

function print_help()
	println("usage:")
	println("statFasta_2.0.jl [-o OUTRGOUP1[,OUTGROUP2,...]] -s <INDIV_STAT_FILE> -f <FASTA_FILES>...\n")
	println("WARNING: <INDIV_STAT_FILE> cannot end with any of the usual fasta file extensions (see 'http://en.wikipedia.org/wiki/FASTA_format#File_extension')")
	println("OUTRGOUPs are race names, not individual names. Regexp are not allowed for the moment")
end

function get_arg(args, i, fun = x->x)
	if length(args) < i+1
		error("expected argument after $(args[i])")
	end
	i += 1
	fun(args[i]), i
end

type Args
	outgroups
	indivstat
	j
end


## II) script
#println("")
#println(length(ARGS))
#println(join(ARGS[1:8], " "))
#println(!ismatch(r"-s ", join(ARGS)))
#println("")

if length(ARGS) < 2 || 
	ismatch(r"\.fasta$|\.fas$|\.fa$|\.seq$|\.fsa$|\.fna$|\.ffn$|\.faa$|\.frn$", ARGS[1]) ||
		!ismatch(r"-s ", join(ARGS, " ")) ||
		!ismatch(r"-f ", join(ARGS, " "))

	println("--------------")
	if !ismatch(r"-s ", join(ARGS, " ")) 
		println("Argument -s is lacking")
	end
	if !ismatch(r"-f ", join(ARGS, " "))
		println("Argument -f is lacking")
	end
	print_help()
	println("--------------")

	exit()
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
println(vecfile, "file\trace\tclone\thaplotype\tNbcalledSites")

# print header of stdout
println("fasta\tlg\tfNpos\tfNposIng\tpNpos\tpNposIng\tfCalled\tfCalledIng\tfCalledClean\tfCalledCleanIng\tfCalledN\tfCalledNIng\tfCalledIns\tfCalledInsIng\tfCalledNIns\tfCalledNInsIng\tfPolym\tfPolymIng\tfPolCle\tfPolCleIng\tfPolN\tfPolNIng\tfPolIns\tfPolInsIng\tfPolNIns\tfPolNInsIng\tpNpos\tpNIng\tfInspos\tfInsposIng\tpInsPos\tpInsPosIng\tpIns\tpInsIng\tfullColN")

for f in args.j:length(ARGS)

	fasta = ARGS[f]

	# apply read_fasta to file, then close it again
	chunks = open(read_fasta, fasta)

	whole, ingrp = compStats(chunks, vecfile, fasta, args.outgroups)
	
	lg, Ncalled, NcalledClean, NcalledN, NcalledIns, NcalledNIns, Npolym, NpolCle, NpolN, NpolIns, NpolNIns, Npos, PcN, fullColN, Inspos, propInsPos, propN, propIns, NcalledIng, NcalledCleanIng, NcalledNIng, NcalledInsIng, NcalledNInsIng, NpolymIng, NpolCleIng, NpolNIng, NpolInsIng, NpolNInsIng, NposIng, PcNIng, InsposIng, propInsPosIng, propNIng, propInsIng = lastStats(whole, ingrp)

	println(basename(fasta), "\t", lg, "\t",
			Npos, "\t", NposIng, "\t",
			PcN, "\t", PcNIng, "\t",
			Ncalled, "\t", NcalledIng, "\t",
			NcalledClean, "\t", NcalledCleanIng, "\t",
			NcalledN, "\t", NcalledNIng, "\t",
			NcalledIns, "\t", NcalledInsIng, "\t",
			NcalledNIns, "\t", NcalledNInsIng, "\t",
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

