#Use this method to identify nullomers from stranded cfRNA. 

using DelimitedFiles, ArgParse, Random, Statistics, BioSequences, StatsBase;

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--nullomer_length", "-l"
        help = "nullomer length"
        arg_type = Int
        default = 12
        "--maxgnomadprob", "-p"
        help = "maximum probability for finding in the germline to allow for if we are going to retain the nullomer"
        arg_type = Float64
        default = 0.05
        "--minqual", "-q"
        help = "minimum quality per base for reads"
        arg_type = Int
        default = 30
        "--readssavefilename", "-S"
        help = "where the reads files are saved, useful for snakemake"
        arg_type = String
        default = nothing
        "--savename", "-s"
        help = "name to use for saving the file, typically corresponds to the patient"
        arg_type = String
        default = ""
        "--file", "-f"
        help = "location of file reads data on disk"
        arg_type = String
        default = ""
        "--tissue", "-t"
        help = "which tissue to look at"
        arg_type = String
        default = ""
        "--patient", "-P"
        help = "which patient to look at"
        arg_type = String
        default = ""
        "--nullomerfile", "-n"
        help = "list of nullomers"
        arg_type = String
        default = ""
        "--logfile"
        help = "file to store metadata"
        arg_type = String
        default = ""
        "--summaryfile"
        help = "file to store summary results"
        arg_type = String
        default = ""
        "--regions", "-R"
        help = "path to .bed-file indicating which regions of the genomes were captured"
        arg_type = String
        default = ""
        "--nullomerpath", "-N"
        help = "path to directory where files are stored"
        arg_type = String
        default = "/lustre/scratch117/cellgen/team218/MH/nullomers/"
        "--proportioncaptured", "-c"
        help = "proportion of reads that are in the capture regions"
        arg_type = Float64
        default = 0.0
    end
    return parse_args(s)
end

function parseRead(nullomerlen, seq, qualscores, minqual, nullomer2ind, nulllistind)
    #Use this method to parse either data from a .bam or .fastq file
    nf = Vector{Any}(undef, nulllistind);
    for nind in 1:nulllistind
        nf[nind]=Int64[];
    end
    lowqualitymatches = Int64[]
    for i in 1:(length(seq)-nullomerlen+1)
                kmer = seq[i:(i+nullomerlen-1)]
        if kmer!=nothing
            for nind in 1:nulllistind
                if haskey(nullomer2ind[nind], kmer)
                    if (qualscores!=nothing) && (minimum(qualscores[i:(i+nullomerlen-1)])>minqual) #check that no poor bases overlap kmer of interest
                        push!(nf[nind], nullomer2ind[nind][kmer])
                    end
                end
            end
        end
    end
    return nf
end


function findNullomersInFastqFile(nullomerpath, nullomerfile, filename::String="", nullomerlen::Int64=16, minqual=13, savereadsname=nothing, genomesize=2945849067)
    #Use this method to find nullomers from a fastq file
    #nullomerindskeep indexes of nullomers to keep
    covtarget = 100000000.0

    nullomerfiles = split(nullomerfile, ",")
    nulllistind=length(nullomerfiles)
    nullomers = Vector{Any}(undef, nulllistind)
    nullomerindskeep = Vector{Any}(undef, nulllistind)
#    nullomerps = Vector{Any}(undef, nulllistind)
    nullomer2ind = Vector{Any}(undef, nulllistind)
    nullsuffix= Vector{Any}(undef, nulllistind)

    for nind in 1:nulllistind
        tmpstrA=split(nullomerfiles[nind],"/");
        nullsuffix[nind]=split(tmpstrA[length(tmpstrA)],"_")[1];
        nullomers[nind] = convert(Array{String, 1}, readdlm(nullomerpath * nullomerfiles[nind], comments=true)[:,1]);
        nullomerindskeep[nind] = 1:length(nullomers[nind])
#        nullomerps[nind] = hcat(nullomerindskeep[nind], zeros(Int, length(nullomers[nind]), 1))
        println("Read " * string(length(nullomers[nind])) * " nullomers from " * nullomerpath * nullomerfiles[nind] * " and keeping " * string(length(nullomerindskeep[nind])))
        nullomer2ind[nind] = Dict{String, Int64}(String(nullomers[nind][i]) => i for i in nullomerindskeep[nind]);
    end

    nullomersfound = Int64[]
    nreads = 0;
    nbps = 0;
    nextlineseq = false;
    nextlinequal = false
    qualline = ""
    headerline = ""
    seq = "";
    savereadsfile = Vector{Any}(undef, nulllistind)
    if savereadsname!=nothing
        for nind in 1:nulllistind
            filenameTmp = savereadsname * "_" * nullsuffix[nind] * ".fastq"
            println("Saving reads with " * nullsuffix[nind] * " nullomers to " * filenameTmp)
            savereadsfile[nind] = open(filenameTmp, "w" )
        end
    end
    cmdstr = `cat $filename`;
    if endswith(filename, ".fastq.gz") || endswith(filename, ".fq.gz")
        cmdstr = `zcat $filename`
    elseif startswith(filename, "SRR") || endswith(filename, ".sra") #warning! should only be done with single-end experiments
        cmdstr = `/home/jovyan/sratoolkit.2.11.0-ubuntu64/bin/fastq-dump -Z $filename`
    end

    for line in eachline(cmdstr) #TODO: Understand this beter
        if length(line)>0
            if (line[1]=='@') && !nextlineseq && !nextlinequal
                nextlineseq = true;
                headerline = line
            elseif nextlineseq
                seq = line
                nextlineseq = false;
            elseif (line[1]=='+') && !nextlineseq && !nextlinequal
                nextlinequal = true;
                qualline = line
            elseif nextlinequal
                qualscores = Int64[];
                if (length(line)>=nullomerlen) && (seq!=nothing) && (length(line)==length(seq))
                    for i in 1:length(line)
                        try
                            push!(qualscores, Int(line[i])-33);
                        catch e
                            println("error with qualscore char")
                        end
                    end
# modify call to take multiple nullomer lists and save reads seperately -------------------------
                    nf = parseRead(nullomerlen, seq, qualscores, minqual, nullomer2ind, nulllistind)
                    for nind in 1:nulllistind
                        if (length(nf[nind])>0)
                            if savereadsfile[nind]!=nothing #also save the reads to a file in case we want to reprocess later. This should greatly speed things up as we do not need to process the entire file.
                                println(savereadsfile[nind] , headerline * "\n" * string(seq) * "\n" * qualline * "\n" * line)
                            end
                            append!(nullomersfound, nf[nind])
                        end
                    end
                    nbps += length(seq)
                end
                nextlinequal = false;
                nreads += 1;
            end
        end
    end
    nullomerinds, nullomercounts = rle(sort(nullomersfound));
    for nind in 1:nulllistind
        if savereadsfile[nind]!=nothing
            close(savereadsfile[nind])
        end
    end
    return nbps, nullomersfound
end

function parseLogFile(filename, readsfile, fieldname)
    patientfound = false;
    value = nothing
    for line in eachline(filename)
        tokens = split(line, ":")
        if patientfound && (strip(tokens[1])==fieldname)
            value = parse(Float64, string(strip(tokens[2])))
            if value>0
                break;
            end
        elseif !patientfound && (strip(tokens[1])=="readsfile")
            patientfound = strip(tokens[2])==readsfile;
        end
    end
    if value==nothing
        println("Warning, could not find coverage for " * readsfile)
        coverage = 0.0;
    else
        println("Found coverage " * string(round(value; digits=2)) * " for " * readsfile)
    end
    return value
end

function findReverseComplement(nullomersstrs, savefilename)
    if typeof(nullomersstrs)==String
        nullomersstrs = convert(Array{String, 1}, readdlm(nullomersstrs, comments=true)[:,1])
    end
    nullomers = LongDNASeq[];
    indsfound = Bool[];
    for n in nullomersstrs
        push!(nullomers, LongDNASeq(n));
        push!(indsfound, false);
    end

    savefilehandle = open(savefilename, "w")
    for i in 1:length(nullomers)
        if !indsfound[i]
            ind = findfirst(nullomersstrs.==string(reverse(complement(nullomers[i]))))
            if ind!=nothing
                println(savefilehandle, string(i) * "\t" * string(ind));
                indsfound[i] = true;
                indsfound[ind] = true;
            end
        end
    end
    close(savefilehandle)
end
function main()
    parsed_args = parse_commandline();
    nullomerlen = parsed_args["nullomer_length"];
    maxbps = Inf
    maxgnomadprob = parsed_args["maxgnomadprob"];
    tissue = parsed_args["tissue"];
    patient = parsed_args["patient"];
    minqual = parsed_args["minqual"]
    savename = parsed_args["savename"]
    readssavefilename = parsed_args["readssavefilename"]
    readsfile = parsed_args["file"]
    nullomerfile = parsed_args["nullomerfile"]
    logfile = parsed_args["logfile"]
    summaryfile = parsed_args["summaryfile"]
    regionsfilename = parsed_args["regions"]
    proportioncaptured = parsed_args["proportioncaptured"]
    nullomerpath = parsed_args["nullomerpath"]

    nullomerliststart = "TargetedMaps/TissueSpecific/"

    genomesize = 2945849067
    hsgenomesize = 2945849067
    regions = nothing
    regionsintervalcollection = nothing
    if regionsfilename!=""
        @time regions = readdlm(nullomerpath * regionsfilename, comments=true);
        regionsintervalcollection = IntervalCollection{}([Interval(regions[i,1], regions[i,2], regions[i,3], '.') for i in 1:size(regions, 1)], true);
        genomesize = sum(regions[:,3] - regions[:,2])
        println("Found " * string(size(regions, 1)) * " regions of length " * string(genomesize) * " bps from " * nullomerpath * regionsfilename * " for an enrichment of " * string(round(hsgenomesize/genomesize; digits=2)))
    end
    captureenrichmentratio = hsgenomesize/genomesize
    fname = nullomerpath * readsfile
#   modified to read multiple nullomer lists  and then thi call below
    nbps, nullomersfound = findNullomersInFastqFile(nullomerpath, nullomerfile, fname, nullomerlen, minqual, readssavefilename, genomesize);
    if readssavefilename!=nothing
        #save the results to a logfile
        logfilehandle = try
            open(nullomerpath * logfile, "a")
        catch
            nothing
        end
        if logfilehandle!=nothing
            println(logfilehandle, "{")
            println(logfilehandle, "\tpatient: " * savename)
            println(logfilehandle, "\treadsfile: " * readsfile)
            println(logfilehandle, "\tcoverage: " * string(round(nbps/genomesize; digits=4)))
            println(logfilehandle, "}")
            close(logfilehandle)
        end
    else
        rcmapfilename = nullomerpath * replace(nullomerfile, ".tsv" => "_rcmap.tsv")
        if (!Base.Filesystem.isfile(rcmapfilename))
            findReverseComplement(nullomerpath * nullomerfile, rcmapfilename);
        end

        rcmap = convert(Array{Int64, 2}, readdlm(rcmapfilename, comments=true));
        rcmapf = Dict{Int64, Int64}();
        rcmapr = Dict{Int64, Int64}();
        for i in 1:size(rcmap, 1)
            if (nullomerps==nothing) || (findfirst(nullomerindskeep.==rcmap[i,1])!=nothing)
                rcmapf[rcmap[i,1]] = rcmap[i,2];
                rcmapr[rcmap[i,2]] = rcmap[i,1];
            end
        end
    #Find all of the substitution mutations that can create this nullomer throughout the genome (or the part that we are considering)
    nullomermapfile = nullomerpath * replace(nullomerfile, ".tsv" => "_subs.tsv")
    nullomermap = Int64[];
    nullomermapcapture = Int64[];
    if Base.Filesystem.isfile(nullomermapfile)
    for line in eachline(nullomermapfile)
        if !startswith(line, "#")
            tokens = split(line, "\t");
            if length(tokens)==5
                overlaps = true
                if regionsintervalcollection!=nothing
                    overlaps = false
                    chrstr = "chr" * string(tokens[1])
                    if tokens[1]==23
                        chrstr = "chrX"
                    elseif tokens[1]==24
                        chrstr = "chrY"
                    end
                    mutloc = Interval(chrstr, parse(Int, tokens[2]), parse(Int, tokens[2])+1, '.')
                    for r in regionsintervalcollection
                        if isoverlapping(r, mutloc)
                            overlaps = true
                            break;
                        end
                    end
                end
                n = tryparse(Int64, tokens[5])
                if findfirst(nullomerindskeep.==n)!=nothing
                    if overlaps & (regionsintervalcollection!=nothing)
                        push!(nullomermapcapture, n);
                    else
                        push!(nullomermap, n);
                    end
                end
            end
        end
    end
    end
    nullomermapind, nullomermapcount = rle(sort(nullomermap));
    nullomermapcaptureind, nullomermapcapturecount = rle(sort(nullomermapcapture));

        genomecov = parseLogFile(nullomerpath * logfile, readsfile, "coverage")
        nullomerinds, nullomercounts = rle(sort(nullomersfound));

        nmutations = vec(zeros(size(nullomerindskeep)));
        for i in 1:length(nullomerindskeep)
            #count the occurrences in the map
            ind = findfirst(nullomerindskeep[i].==nullomermapind);
            if (ind!=nothing) && (length(nullomermapcount)>=ind) && (nullomermapcount[ind]!=nothing)
                nmutations[i] = (1-proportioncaptured)*genomecov*nullomermapcount[ind];
            end
            ind2 = nothing;
            if haskey(rcmapf, nullomerindskeep[i])
                ind2 = findfirst(rcmapf[nullomerindskeep[i]].==nullomermapind);
            elseif haskey(rcmapr, nullomerindskeep[i])
                ind2 = findfirst(rcmapr[nullomerindskeep[i]].==nullomermapind);
            end
            if ind2!=nothing
                nmutations[i] += (1-proportioncaptured)*genomecov*nullomermapcount[ind2];
            end
            #now count the mutations that were found in the capture regions. These get inflated since the number of
            ind = findfirst(nullomerindskeep[i].==nullomermapcaptureind);
            if ind!=nothing
                nmutations[i] += genomecov*captureenrichmentratio*proportioncaptured;
            end
            ind2 = nothing;
            if haskey(rcmapf, nullomerindskeep[i])
                ind2 = findfirst(rcmapf[nullomerindskeep[i]].==nullomermapcaptureind);
            elseif haskey(rcmapr, nullomerindskeep[i])
                ind2 = findfirst(rcmapr[nullomerindskeep[i]].==nullomermapcaptureind);
            end
            if ind2!=nothing
                nmutations[i] += genomecov*captureenrichmentratio*proportioncaptured;
            end
        end
         #lambdas, bginds, somaticinds, nnullomers, lambdasgerm = analyzeGenomicNullomers(nullomerps, nullomerindskeep, nullomerinds, nullomercounts, maxallelefreq, nmutations, genomecov, fdr, pairedend; nullomerlen=nullomerlen, errorrate=errorrate);
        nnullomers = Int[]
        n = length(nullomerindskeep)
        for i in 1:n
            #count the occurrences in the data, rcs have already been included by findNullomersInFastqFile
            ind = findfirst(nullomerindskeep[i].==nullomerinds);
            if ind!=nothing
                push!(nnullomers, nullomercounts[ind]);
            else
                push!(nnullomers, 0);
            end
        end
        #save all of the nullomers that were found into a file
        filehandle = open(nullomerpath * savename, "w")
        println(filehandle, "#neomer\tnreads\tlambda\tnmutations\tp_germline\tallelefreq")
            for i in 1:n
                if nnullomers[i]>0
                    println(filehandle, string(nullomerindskeep[i]) * "\t" * string(nnullomers[i]) * "\t" * string(round(lambda[i]; digits=5)) * "\t" * string(round(nmutations[i]/genomecov; digits=3)) * "\t" * string(round(nullomerps[nullomerindskeep[i],2]; digits=6)) * "\t" * string(round(a[i]; digits=5)))
                end
            end
        close(filehandle)
        #save a summary that can be used for plotting
        if !isfile(nullomerpath * summaryfile) #write a header
            filehandle = open(nullomerpath * summaryfile, "w")
            println(filehandle, "#sample\ttissue\tcoverage\tnullomerlength\tmaxgnomadprob\terrorrate\tmax_allele_freq\tprop_captured\ttotal_no_nullomers\tno_nullomer_reads\tallelefreq")
            close(filehandle)
        end
        filehandle = open(nullomerpath * summaryfile, "a")
        println(filehandle, patient * "\t" * tissue * "\t" * string(round(genomecov; digits=2)) * "\t" * string(nullomerlen) * "\t" * string(round(maxgnomadprob; digits=5)) * "\t" * string(round(errorrate; digits=4)) * "\t" * string(round(maxallelefreq; digits=4)) * "\t" * string(round(proportioncaptured; digits=4)) * "\t" * string(length(nullomerindskeep)) * "\t" * string(sum(nnullomers)) * "\t" * string(round(sum(a)/n; digits=5)));
        close(filehandle)
    end
end

main()
