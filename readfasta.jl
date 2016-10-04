
type Chunk
	id
	scaffold
	pos_min
	pos_max
	strand
	data
end

# constructor functions have the same name as types
function Chunk(line)
	# split at # and remove whitespace
	fields = map(strip, split(line[2:end], '#'))
	# julia starts counting at 1
	id = fields[1]
	# get scaffold and coordinates
	sc_pos = split(fields[2], ':')
	scaffold = sc_pos[1]
	poss = map(int, split(sc_pos[2], '-'))
	
	# call built-in constructor
	Chunk(id, scaffold, poss[1], poss[2], fields[3][end], Char[])
end

# add a line
function add_data(chunk, line)
	# collect converts string into array of char
	append!(chunk.data, collect(line))
end

# get a single letter
function get_at_pos(chunk, pos::Int)
	if pos < chunk.pos_min || pos > chunk.pos_max
		error("trying to access position $pos in chunk $(chunk.scaffold):$(chunk.pos_min)-$(chunk.pos_max)")
	end

	chunk.data[pos-chunk.pos_min + 1]
end

# get a range of letters (as array of char)
function get_at_pos(chunk, range::Range)
	if first(range) < chunk.pos_min || last(range) > chunk.pos_max
		println("WARNING: trying to access range $(first(range)):$(last(range)) " * 
			"in chunk $(chunk.scaffold):$(chunk.pos_min)-$(chunk.pos_max)")
	end

	# adjust requested range to available positions
	p1 = max(first(range)-(chunk.pos_min-1), 1)
	p2 = min(last(range)-(chunk.pos_min-1), length(chunk.data))

	chunk.data[p1:p2]
end

# find the chunk that belongs to a position
function find_at_pos(chunks, pos)
	# we don't have that one
	if pos < chunks[1].pos_min || pos > chunks[end].pos_max
		return 0
	end

	# chunks are sorted by first pos, so just do a simple binary search

	min = 1
	max = length(chunks)

	while max >= min
		test = div(min+max, 2)

		if pos < chunks[test].pos_min
			max = test - 1
		elseif pos > chunks[test].pos_max
			min = test + 1
		else
			return test
		end
	end

	return 0
end


function read_fasta(input)
	# scaffold => chunk
	chunks = Dict{ASCIIString, Vector{Chunk}}()

	scaffold = "" 

	for line in eachline(input)
		# strip removes whitespace from beginning and end
		line = strip(line)

		# header
		if line[1] == '>'
			# create a chunk from the line
			chunk = Chunk(line)
			scaffold = chunk.scaffold

			# if the scaffold is new we have to add a list to the dict
			if ! haskey(chunks, scaffold)
				chunks[scaffold] = Chunk[]
			end
			# add chunk to end of list for this scaffold
			push!(chunks[scaffold], chunk)

		# data
		else
			add_data(chunks[scaffold][end], line)
		end

		if eof(input)
			# add the last one in the file
			push!(chunks[scaffold], chunk)
		end
	end


	# sort chunks by start position
	for chunk_list in values(chunks)
		sort!(chunk_list, by = c -> c.pos_min)
	end
	
	chunks
end

