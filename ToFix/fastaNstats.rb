#!/usr/bin/env ruby

# count number of N in sequence from stdin
# does not work right now

class Chunk
	attr_reader :data

	def initialize
		@data = ""
	end

	def add_data(line)
		@data += line
	end
end

def read_fasta(input)
	chunks = []

	input.readlines.each do |line|
		line.strip!

		if line.empty?
			next
		end

		if line.starts_with?(">")
			chunks << Chunk.new
		else
			chunks.last.add_data(line)
		end
	end
end

def stats(chunks)
	hasN = [false] * chunks[0].data.length

	countN = 0

	chunks.each do |chunk|
		if chunk.data.length != hasN.length
			raise 'inconsistent array lengths!'
		end

		for i in 0..chunk.data.length-1 do
			if chunk.data[i] == 'N'
				countN += 1
				hasn[i] = true
			end
		end
	end
    
    prop = countN/(hasN.length*chunks.length)

	hasN, prop
end

chunks = read_fasta($stdin)

hasN, propN = stats(chunks)

puts(propN)

puts(hasN.join(", "))
