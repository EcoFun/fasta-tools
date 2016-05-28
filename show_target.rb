#!/usr/bin/env ruby
# usage:
    # cat my_fasta_file.fasta | show_target.rb REGEX_locus-name > results


match = Regexp.new(ARGV[0])

do_print = false
first = true

$stdin.each_line do |line|
	if line.start_with?(">")
		do_print = false
		if line =~ match
			do_print = true
            if first == true
                puts line
            end
        first = false
		end
	elsif do_print
		print line.gsub("\n","")
	end
end
