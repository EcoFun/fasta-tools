#!/usr/bin/env julia

# dependencies
# - julia version 0.4.5

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

usage="
SYNOPSIS:
split_fasta.jl [-h] [-b] [-D] [-r] [-v] [-d output_dir] [-x pattern] [-y string] fasta_file

DESCRIPTION:
Split a multi-fasta file into as many files as sequences contained.
The names of the output files will be the headers. Thus, some options 
allow to replace problematic characters (see OPTIONS).

WARNING: the short options cannot be concatenated together! (e.g. '-Ddr' is invalid, see example below)

ARGUMENTS:
    fasta_file      input fasta file

OPTIONS:
    -b              do not replace blanck by underscore in headers (done by default)
    -d              output directory
    -D              debug mode
    -h              print this help and exit (irrespectively to other arguments)
    -r              allow regular expression for pattern replacing? (default is false)
                    Works with option '-x'
    -v              verbose mode
    -x              pattern to be removed in headers (future file name)
    -y              string to replace '-x' patterns with

EXAMPLE:
split_fasta.jl -r -v './2016-2017_Venturia/Data4Corentin/00_Data/Single_CDS' -p mRNAVein -r B04 './2016-2017_Venturia/Data4Corentin/00_Data/B04/cds'
"

# functions
function get_arg(args, i, fun = x->x)
	if length(args) < i+1
		error("expected argument after $(args[i])")
	end
	i += 1
	fun(args[i]), i
end

# get parameters
i = 1 ; 
dir = "." ; pat1 = "" ; pat2 = "" ; inp = ""
sprep = "true" ; debug = "false" ; reg= "false" ; verbose= "false"
gooda = r"-[bdDhrvxy]"
while i ≤ length(ARGS)
    a = ARGS[i]
    
    if i != length(ARGS) && !ismatch(gooda, a)
        println("ERROR: Bad argument '$a' in the command line")
        println("       !!! Beware that concatenation of options is not allowed (e.g. '-Ddr' is invalid) !!!")
        println(usage)
        exit()
    end
    
#~    println(a)
    if a == "-h"          # print help then exit
        println(usage)
        exit()
    elseif a == "-b"    # replace blanks by underscore in header?
        sprep = "false"
    elseif a == "-d"    # directory where to write the outputs
        dir, i = get_arg(ARGS, i)
        if ismatch(gooda, dir)
            println("ERROR: bad value '$dir' for flag '$a'")
            println(usage)
            exit()
        end
    elseif a == "-r"    # allow regular expression?
        reg = "true"
    elseif a == "-x"    # pattern to be removed in seq name
        pat1, i = get_arg(ARGS, i)
        if ismatch(gooda, pat1)
            println("ERROR: bad value '$pat1' for flag '$a'")
            println(usage)
            exit()
        end
    elseif a == "-y"    # pattern to replace '-p' with
        pat2, i = get_arg(ARGS, i)
        if ismatch(gooda, pat2)
            println("ERROR: bad value '$pat2' for flag '$a'")
            println(usage)
            exit()
        end
    elseif a == "-D"    # debbug
        debug = "true"
    elseif a == "-v"    # verbose mode
        verbose = "true"
    else
        inp = ARGS[end] # input file
    end
    i += 1
end

if debug == "true"
    println(ARGS)
    println("-b=$sprep; -d=$dir; -D=$debug; -r=$reg; -x=$pat1; -y=$pat2; -v=$verbose; input=$inp")
    exit()
end


## real thing
# set variables and open input file
println("Process file '$inp'")
run(`mkdir -p $dir`)
f = open(inp)

i = 0
nom =""
for l in eachline(f)
    if l[1] == '>'
        # close previous fasta file if needed
        if i != 0
            close(fas)
        else
            i = 1
        end
        
        # set new fasta file name and open it
        nom = strip(replace(l, r"[>+:-]", ""))
        if sprep == "true"
            nom = replace(nom, " ", "_")
            nom = replace(nom, "\t", "__")
        end
        if length(pat1) > 0
            if reg == "true"
                nom = replace(nom, Regex(pat1), pat2)
            else
                pat1 = replace(pat1, r"[>+:-]", "")
                nom = replace(nom, pat1, pat2)
            end
        end
        out = dir * "/" * nom * ".fas"
        global fas = open(out, "w")

        # print line
        if verbose == "true"
            println("Write '$l' in file '$out'")
        end
        print(fas, l)
    else
        print(fas, strip(l))
    end
end

close(f)
