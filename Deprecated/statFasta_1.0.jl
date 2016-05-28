#!/usr/bin/env julia0.3.11

#### script analysing the information available in the multifasta files
# This version is only able to compute the stats for the whole dataset but is very quick

#!/usr/bin/env julia

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


# calculate stats
function stats(chunks::Vector{Chunk}, vecfile, fasta)

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
	
	ref = []
	
	for c in 1:length(chunks)
		chunk = chunks[c]

		# make sure all chunks are the same length
		@assert	length(chunk.data) == length(hasN)
		
		if c == 1
			ref = chunk.data
		end
		
		# start incrementation individual informative sites
		indi_info = 0
		
		# 1) check all sites
		for i in 1:length(chunk.data)

			# 1.1) count number of Ns
			if chunk.data[i] == 'N'
				countN += 1
				hasN[i] = 1

			# 1.2) count number of indel positions
			elseif chunk.data[i] == '-'
				countIns += 1
				hasIns[i] = 1

			# 1.3) record informative sites
			elseif chunk.data[i] == 'A' || chunk.data[i] == 'T' || 
				chunk.data[i] == 'G' || chunk.data[i] == 'C'
				
				# 1.3.1) count number of informative sites for individual i
				indi_info += 1

				# 1.3.2) count number of polymorphic sites
					# 1.3.2.1) update ref sequence
				if (ref[i] == 'N' || ref[i] == '-')
					ref[i] = chunk.data[i]
				end

					# 1.3.2.2) check polymorphism
				if chunk.data[i] != ref[i]
					polym[i] = 1
				end
			else
				error("Invalid symbol" * chunk.data[i] * " in alignment: file "
					 * fasta * ", haplotype " * c * " (valid symbols are 'A', 'C', 'G', 'T', 'N' and '-')")
			end

			if i == length(chunk.data)
				println(vecfile, fasta, "\t", chunk.race, "\t", chunk.nom, "\t",  chunk.haplo,"\t",  indi_info) 
				indi_info = 0
			end

			# 1.4) if last chunk, perform few more tests
			if c == length(chunks)
				# 1.4.0) check whether the site is informative (at least one individual with no N nor indel)
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
				
				# 1.4.1) check whether the site is clean (no Ns, no indels) - if so necessarily informative
				if hasN[i] == 0 && hasIns[i] == 0
					infoClean[i] = 1
				end
				
				# 1.4.2) for polymorphic sites:
				if polym[i] == 1
					# 1.4.2.1) count the number of clean polym
					if hasN[i] == 0 && hasIns[i] == 0
						polCle[i] = 1
					end
					
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

			end	# end of 'if c == length(chunks)'
		end	# end of 'for i in 1:length(chunk.data)'
	end	# end of 'for c in 1:length(chunks)'

	# functions can return several values
	propN = countN/(length(hasN)*length(chunks))
	propIns = countIns/(length(hasIns)*length(chunks))
	
	hasN, propN, hasIns, propIns, infoClean, polym, polCle, polN, polIns, polNIns, info, infoN, infoIns, infoNIns
end

function print_help()
	println("--------------")
	println("usage:")
	println("fastaNstats.jl <INDIV_STAT_FILE> <FASTA_FILES>...\n")
	println("WARNING: <INDIV_STAT_FILE> cannot cannot end with any of the usual fasta file extensions (see 'http://en.wikipedia.org/wiki/FASTA_format#File_extension')")
	println("--------------")
end


## II) script
if length(ARGS) < 2 || ismatch(r"\.fasta$|\.fas$|\.fa$|\.seq$|\.fsa$|\.fna$|\.ffn$|\.faa$|\.frn$", ARGS[1])
	print_help()
	exit()
end

# open INDIV_STAT_FILE file
vecfile = open(ARGS[1], "w")
println(vecfile, "file\trace\tclone\thaplotype\tNbInfoSites")

# print header of stdout
println("file\tlength\tNbInformative\tNbInfoClean\tNbInfoN\tnbInfoInd\tNbInfoNInd\tNbPol\tNbPolClean\tNbPolN\tNbPolInd\tNbPolNInd\tNbNsites\tPropNsites\tPropNs\tfullColN\tNbIndSites\tPropIndSites\tPropInd")

for f in 2:length(ARGS)

	fasta = ARGS[f]

	# apply read_fasta to file, then close it again
	chunks = open(read_fasta, fasta)

	hasN, propN, hasIns, propIns, infoClean, polym, polCle, polN, polIns, polNIns, info, infoN, infoIns, infoNIns = stats(chunks, vecfile, fasta)

	lg = length(hasN)
	
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
	propNpos = Npos/lg
	fullColN = propN==propNpos

	Inspos = sum(hasIns)
	propInsPos = Inspos/lg

	println(fasta, "\t", lg, "\t", Ninfo, "\t", NinfoClean, "\t", NinfoN, "\t", NinfoIns, "\t", NinfoNIns, "\t", Npolym, "\t", NpolCle, "\t", NpolN, "\t", NpolIns, "\t", NpolNIns, "\t", Npos, "\t", propNpos, "\t", propN, "\t", fullColN, "\t", Inspos, "\t", propInsPos, "\t", propIns)
end

# close outfile
close(vecfile)

