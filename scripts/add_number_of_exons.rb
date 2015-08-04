usage= <<-eof
USAGE: TODO
  ruby add_number_of_exons.rb file.gtf in.csv out.csv
  ____________________________________________________
  Adds number of exons to in.csv for ENMUSG
  IN: file.gtf in.csv
  OUT: out.csv
  in.csv needs to look like this:
    ensgene,pExon,pGroup,pInteraction,fdr
    ENSMUSG00000027365,0,7.72E-11,1.34E-76,2.49E-72
    ENSMUSG00000045103,0,3.55E-05,1.48E-70,2.75E-66
    ENSMUSG00000039219,0,1.43E-13,1.80E-64,3.34E-60
    ENSMUSG00000020741,0,2.92E-41,3.63E-61,6.75E-57
    ENSMUSG00000004056,0,1.02E-24,1.08E-57,2.01E-53
    ENSMUSG00000009376,4.94E-187,2.76E-10,7.70E-53,1.43E-48
  ____________________________________________________

eof

puts usage if ARGV.length != 3

num_of_exons = {}
exons = {}
File.open(ARGV[0]).each do |line|
  line.chomp!
  line = line.split("\t")
  next unless line[2] == "exon"
  line[-1] =~ /(ENSMUSG\d*)/
  name = $1
  exons[name] ||= []
  exons[name] << "#{line[0]}:#{line[3]}-#{line[4]}" unless exons[name].include?("#{line[0]}:#{line[3]}-#{line[4]}")
  #num_of_exons[name] ||= 0
  #num_of_exons[name] += 1
end

exons.each_pair do |key,value|
  num_of_exons[key] = value.length
end

#puts num_of_exons["ENSMUSG00000025909"] # should be 124

out_file = File.open(ARGV[2],'w')
first = true
File.open(ARGV[1]).each do |line|
  line.chomp!
  #puts line
  line = line.split(",")
  if first
    line << "num_of_exons"
    first = false
  else
    line << num_of_exons[line[0]]
  end
  out_file.puts line.join(",")
end

out_file.close()
